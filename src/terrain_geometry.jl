"""Native-free bounded terrain geometry with immutable axes and fields.

Use `TerrainGeometry.geometry` for queries within an explicitly padded served
domain. Geometry alone does not establish atmospheric validity.
"""
module TerrainGeometry
export Ellipsoid, ServedDomain, GeometryTile, geometry, radial_position, padding_certificate

struct Ellipsoid
    equatorial_m::Float64
    polar_m::Float64
    function Ellipsoid(a::Real,b::Real)
        x,y=Float64(a),Float64(b)
        isfinite(x) && isfinite(y) && 0<y<=x || throw(ArgumentError("Finite positive oblate ellipsoid radii required"))
        new(x,y)
    end
end

struct ServedDomain
    geodetic_latitude_deg::Tuple{Float64,Float64}
    east_longitude_deg::Tuple{Float64,Float64}
    ellipsoid_height_m::Tuple{Float64,Float64}
    function ServedDomain(latitude,longitude,height)
        pairs=map((latitude,longitude,height)) do values
            length(values)==2 || throw(ArgumentError("Two domain endpoints required"))
            lo,hi=Float64.(values)
            isfinite(lo) && isfinite(hi) && lo<hi || throw(ArgumentError("Finite increasing domain bounds required"))
            (lo,hi)
        end
        -89<=pairs[1][1]<pairs[1][2]<=89 || throw(ArgumentError("Served latitude must stay within [-89, 89] degrees"))
        pairs[3][1]>=-100000 || throw(ArgumentError("Served height must be at least -100000 metres"))
        0<=pairs[2][1]<pairs[2][2]<360 || throw(ArgumentError("Seam-crossing domains are not supported"))
        new(pairs...)
    end
end

function radial_position(e::Ellipsoid,latitude_deg::Real,height_m::Real)
    p,h=Float64(latitude_deg),Float64(height_m)
    isfinite(p) && -90<=p<=90 && isfinite(h) || throw(ArgumentError("Finite valid geodetic position required"))
    h>e.polar_m^2/e.equatorial_m*-1 || throw(DomainError(h,"Height outside monotone normal-coordinate chart"))
    theta=deg2rad(p);s,c=sin(theta),cos(theta)
    e2=1-(e.polar_m/e.equatorial_m)^2
    n=e.equatorial_m/sqrt(1-e2*s*s)
    x=(n+h)*c;z=(n*(1-e2)+h)*s
    (radius_m=hypot(x,z),planetocentric_latitude_deg=rad2deg(atan(z,x)))
end

# Float64 leaves and tuple storage make the complete object graph immutable.
# Lengths are kept out of the field types to avoid specializing GeometryTile on
# every grid shape. Returned copies are ordinary, independent Julia arrays.
struct _FrozenVector <: AbstractVector{Float64}
    data::Tuple{Vararg{Float64}}
end
Base.size(v::_FrozenVector) = (length(v.data),)
Base.IndexStyle(::Type{_FrozenVector}) = IndexLinear()
Base.getindex(v::_FrozenVector, i::Int) = v.data[i]
Base.copy(v::_FrozenVector) = collect(v.data)

struct _FrozenMatrix <: AbstractMatrix{Float64}
    data::Tuple{Vararg{Float64}}
    dims::Tuple{Int,Int}
end
Base.size(m::_FrozenMatrix) = m.dims
Base.IndexStyle(::Type{_FrozenMatrix}) = IndexLinear()
Base.getindex(m::_FrozenMatrix, i::Int) = m.data[i]
Base.copy(m::_FrozenMatrix) = reshape(collect(m.data), m.dims)

"""
    GeometryTile(latitude, longitude, areoid, terrain, ellipsoid, domain;
                 minimum_padding_deg)

Construct a bounded native-free geometry tile from strictly increasing
planetocentric latitude and east-longitude axes in degrees, and latitude-major
matrices of areoid radius and terrain height in metres. Storage is immutable;
inputs are snapshotted and `copy(tile.areoid_radius_m)` returns an editable copy.
The served domain uses geodetic latitude and ellipsoid height. Declared padding
is required in each axis, separately from the numerical guard. Poles and seam
crossing are unsupported. This is a geometry contract, not atmosphere validity.
"""
struct GeometryTile
    latitude_deg::_FrozenVector
    longitude_deg::_FrozenVector
    areoid_radius_m::_FrozenMatrix
    terrain_height_m::_FrozenMatrix
    ellipsoid::Ellipsoid
    domain::ServedDomain
    required_padding_deg::NTuple{2,Float64}
    actual_padding_deg::NTuple{4,Float64}
    numerical_floor_deg::NTuple{2,Float64}
    guarded_padding_deg::NTuple{4,Float64}
    function GeometryTile(latitude,longitude,areoid,terrain,ellipsoid::Ellipsoid,domain::ServedDomain; minimum_padding_deg)
        lat,lon=Float64.(collect(latitude)),Float64.(collect(longitude))
        for axis in (lat,lon)
            length(axis)>=2 && all(isfinite,axis) && all(diff(axis).>0) || throw(ArgumentError("Finite strictly increasing axes required"))
        end
        -90<first(lat)<last(lat)<90 || throw(ArgumentError("Polar tiles unsupported"))
        0<=first(lon)<last(lon)<360 || throw(ArgumentError("Seam-crossing tiles unsupported"))
        a,t=Matrix{Float64}(areoid),Matrix{Float64}(terrain)
        size(a)==size(t)==(length(lat),length(lon)) || throw(DimensionMismatch("Latitude-major field shape required"))
        all(isfinite,a) && all(isfinite,t) && all(a.>0) || throw(ArgumentError("Finite fields and positive areoid radii required"))
        # In this chart, radial latitude increases with geodetic latitude;
        # at fixed latitude it is monotone in height. Rectangle extrema are corners.
        extrema=[radial_position(ellipsoid,p,h).planetocentric_latitude_deg for p in domain.geodetic_latitude_deg for h in domain.ellipsoid_height_m]
        first(lat)<minimum(extrema)<=maximum(extrema)<last(lat) || throw(ArgumentError("Served domain must lie strictly inside padded latitude axis"))
        first(lon)<domain.east_longitude_deg[1]<domain.east_longitude_deg[2]<last(lon) || throw(ArgumentError("Served domain must lie strictly inside padded longitude axis"))
        length(minimum_padding_deg)==2 || throw(ArgumentError("Declare latitude and longitude minimum padding in degrees"))
        required=Tuple(Float64.(minimum_padding_deg))
        all(x->isfinite(x)&&x>0,required) || throw(ArgumentError("Finite positive minimum padding is required"))
        # This explicit policy rejects near-ulp padding. It is not a proved global
        # error bound for every ellipsoid/chart, nor a physical terrain margin.
        floors=(64eps(max(1.,maximum(abs,lat),maximum(abs,extrema))),
                64eps(max(1.,maximum(abs,lon),maximum(abs,domain.east_longitude_deg))))
        all(required[i]>=floors[i] for i in 1:2) || throw(ArgumentError("Declared padding is below the numerical guard floor"))
        actual=(minimum(extrema)-first(lat),last(lat)-maximum(extrema),
                domain.east_longitude_deg[1]-first(lon),last(lon)-domain.east_longitude_deg[2])
        # Reserve the floor separately from the requested usable padding.
        # Round the subtraction down one representable step. The floor remains
        # an empirical numerical policy, not a universal libm error bound.
        guarded=ntuple(i->prevfloat(actual[i]-floors[i<=2 ? 1 : 2]),4)
        all(guarded[i]>=required[i<=2 ? 1 : 2] for i in 1:4) || throw(ArgumentError("Stored tile does not meet the declared minimum padding after numerical guard"))
        new(_FrozenVector(Tuple(lat)),_FrozenVector(Tuple(lon)),
            _FrozenMatrix(Tuple(vec(a)),size(a)),_FrozenMatrix(Tuple(vec(t)),size(t)),
            ellipsoid,domain,required,actual,floors,guarded)
    end
end

"""Report construction-time angular padding for an immutable tile.

Actual and guarded values are south, north, west, east; requirements/floors are
latitude, longitude. Acceptance uses guarded values. This is a numerical geometry
policy, not atmospheric validity or a proved libm-independent error enclosure.
"""
padding_certificate(tile::GeometryTile)=(required_deg=tile.required_padding_deg,
    actual_south_north_west_east_deg=tile.actual_padding_deg,
    numerical_floor_deg=tile.numerical_floor_deg,
    guarded_south_north_west_east_deg=tile.guarded_padding_deg)

function bracket(axis,q)
    isfinite(q) && first(axis)<=q<=last(axis) || throw(DomainError(q,"Query outside stored geometry tile"))
    i=min(searchsortedlast(axis,q),length(axis)-1)
    i,(q-axis[i])/(axis[i+1]-axis[i])
end
function radial_surface(tile::GeometryTile,latitude_deg::Real,east_deg::Real)
    i,u=bracket(tile.latitude_deg,Float64(latitude_deg));j,v=bracket(tile.longitude_deg,Float64(east_deg))
    weights=((i,j,(1-u)*(1-v)),(i+1,j,u*(1-v)),(i,j+1,(1-u)*v),(i+1,j+1,u*v))
    interpolate(field)=sum(w*field[r,c] for (r,c,w) in weights if w>0)
    (areoid_radius_m=interpolate(tile.areoid_radius_m),terrain_height_m=interpolate(tile.terrain_height_m))
end

"""Return radial geometry in metres for a geodetic query inside the served domain.

Longitude is east-positive and bounded, never implicitly wrapped. Heights are above
an oblate ellipsoid. Clearance is radius minus areoid radius minus terrain height;
it is neither geodetic terrain altitude nor an atmosphere-validity guarantee.
"""
function geometry(tile::GeometryTile,latitude_deg::Real,east_deg::Real,height_m::Real)
    p,l,h=Float64(latitude_deg),Float64(east_deg),Float64(height_m)
    for (q,bounds) in zip((p,l,h),(tile.domain.geodetic_latitude_deg,tile.domain.east_longitude_deg,tile.domain.ellipsoid_height_m))
        isfinite(q) && bounds[1]<=q<=bounds[2] || throw(DomainError(q,"Query outside served geometry domain"))
    end
    position=radial_position(tile.ellipsoid,p,h)
    surface=radial_surface(tile,position.planetocentric_latitude_deg,l)
    (;position...,surface...,areoid_height_m=position.radius_m-surface.areoid_radius_m,
      radial_clearance_m=position.radius_m-surface.areoid_radius_m-surface.terrain_height_m)
end
end
