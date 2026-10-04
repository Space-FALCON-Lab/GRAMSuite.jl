"""Native-free, immutable endpoint winds for a bounded frozen Mars atmosphere.

This component supplies no data and does not select a default atmosphere.
Its current coordinate envelope is the independently checked local Mars case.
"""
module MarsWindEndpoints
using LinearAlgebra
using ..TerrainGeometry
export FieldID, Observation, HarmonicChart, EndpointField, WindSnapshot,
       cardinal_weights, field_weights, required_fields, endpoint, winds

const CONDITION_CAP = 1000.0
const AMPLIFICATION_CAP = 10.0
const PLACEMENT_TOLERANCE_M = 1e-5

struct FieldID
    datum::Symbol
    level_m::Int
    piece::Int
    function FieldID(datum::Symbol, level_m::Integer, piece::Integer=0)
        valid = datum === :areoid ? (level_m % 1000 == 0 && piece == 0) :
                datum === :clearance ? (level_m in (5,30) && piece in (1,2)) : false
        valid || throw(ArgumentError("Use areoid integer-km fields or 5/30 m clearance pieces 1/2"))
        new(datum,Int(level_m),Int(piece))
    end
end

"""One physical endpoint observation. Geometry and wind values use metres and m/s.
Coordinates are planetocentric latitude and east longitude in degrees. Exposure
is checked from the supplied observation geometry, whose provenance is the
snapshot provider's responsibility.
"""
struct Observation
    field::FieldID
    latitude_deg::Float64
    longitude_deg::Float64
    areoid_height_m::Float64
    terrain_height_m::Float64
    horizontal_m_s::NTuple{2,Float64}
    function Observation(field::FieldID,p,l,h,s,wind)
        length(wind)==2 || throw(ArgumentError("Two horizontal wind components required"))
        p,l,h,s=Float64.((p,l,h,s));uv=Tuple(Float64.(wind))
        all(isfinite,(p,l,h,s,uv...)) && -90<p<90 && 0<=l<360 || throw(ArgumentError("Finite observation and valid coordinates required"))
        firstlevel=1000floor(s/1000+1.3)
        exposed=field.datum===:areoid ? (firstlevel<=field.level_m && h-s>30 && abs(h-field.level_m)<=PLACEMENT_TOLERANCE_M) : abs(h-s-field.level_m)<=PLACEMENT_TOLERANCE_M
        exposed || throw(ArgumentError("Observation does not expose the named physical field"))
        new(field,p,l,h,s,uv)
    end
end

function series_coefficients(n,s)
    polys=[Float64[1],Float64[0,1]]
    for k in 2:n-1
        p=zeros(k+1)
        for (i,v) in enumerate(polys[end]);p[i+1]+=2v;end
        for (i,v) in enumerate(polys[end-1]);p[i]-=v;end
        push!(polys,p)
    end
    terms=n==5 ? ((3,5s*s),(1,4s^4)) : ((8,10s*s),(6,33s^4),(4,40s^6),(2,16s^8))
    Tuple(map(polys) do p
        c=vcat(p,zeros(48-length(p)))
        for j in 0:47-n
            c[j+n+1]=-sum(v*c[j+k+1]/prod(j+k+1:j+n) for (k,v) in terms)
        end
        Tuple(c)
    end)
end
function polynomial(c,x)
    y=0.0
    for v in Iterators.reverse(c);y=y*x+v;end
    y
end
function basis_values(coeff,lo,hi,l)
    t=(2Float64(l)-lo-hi)/(hi-lo)
    isfinite(t) && abs(t)<=2 || throw(DomainError(l,"Outside chart extrapolation bound"))
    map(c->polynomial(c,t),coeff)
end

"""Well-conditioned scaled harmonic chart; stored arrays are immutable tuples.
Five functions represent lower fields, ten represent surface pieces. Optional
explicit Chebyshev nodes preserve the generator's recorded Float64 coordinates.
"""
struct HarmonicChart
    lo::Float64
    hi::Float64
    size::Int
    nodes::Tuple{Vararg{Float64}}
    coefficients::Tuple
    inverse::Tuple
    condition::Float64
    function HarmonicChart(lo,hi,n::Integer;nodes=nothing)
        lo,hi=Float64.((lo,hi))
        isfinite(lo) && isfinite(hi) && 0<=lo<hi<360 && hi-lo<=9 && n in (5,10) || throw(ArgumentError("Finite bounded 5/10-function chart required"))
        expected=ntuple(i->lo+(hi-lo)*(1-cos(pi*(i-1)/(n-1)))/2,n)
        points=nodes===nothing ? expected : Tuple(Float64.(nodes))
        length(points)==n && all(isfinite,points) && all(points[i]<points[i+1] for i in 1:n-1) && all(abs(points[i]-expected[i])<=1e-12 for i in 1:n) || throw(ArgumentError("Recorded increasing Chebyshev nodes required"))
        c=series_coefficients(n,deg2rad((hi-lo)/2))
        a=reduce(vcat,[permutedims(collect(basis_values(c,lo,hi,l))) for l in points])
        inverse=inv(a);condition=opnorm(a,Inf)*opnorm(inverse,Inf)
        isfinite(condition) && condition<=CONDITION_CAP || throw(ArgumentError("Chart condition cap exceeded"))
        new(lo,hi,n,points,c,Tuple(Tuple(row) for row in eachrow(inverse)),condition)
    end
end
function cardinal_weights(c::HarmonicChart,l::Real)
    b=basis_values(c.coefficients,c.lo,c.hi,l)
    w=ntuple(i->sum(b[j]*c.inverse[j][i] for j in 1:c.size),c.size)
    abs(sum(w)-1)<=1e-10 && sum(abs,w)<=AMPLIFICATION_CAP || throw(DomainError(l,"Chart weight amplification cap exceeded"))
    w
end

struct EndpointField
    id::FieldID
    latitudes::NTuple{2,Float64}
    charts::NTuple{2,HarmonicChart}
    observations::Tuple{Vararg{Observation}}
    function EndpointField(id::FieldID,lats,charts,observations)
        length(lats)==length(charts)==2 || throw(ArgumentError("Two fitting rows required"))
        lat=Tuple(Float64.(lats));cs=Tuple(charts);obs=Tuple(observations)
        all(isfinite,lat) && 33.75<=lat[1]<lat[2]<=37.5 || throw(ArgumentError("Fitting rows must share the supported Mars latitude cells"))
        n=id.datum===:areoid ? 5 : 10
        all(c->c isa HarmonicChart && c.size==n,cs) && length(obs)==2n || throw(ArgumentError("Missing or incompatible endpoint observations"))
        for (row,c) in enumerate(cs)
            bounds=id.datum===:areoid ? (40.0,42.5) : id.piece==1 ? (40.0,40.5) : (40.5,42.5)
            bounds[1]<=c.lo<c.hi<=bounds[2] || throw(ArgumentError("Chart crosses the declared surface piece or supported region"))
            for i in 1:n
                o=obs[(row-1)*n+i]
                o isa Observation && o.field==id && o.latitude_deg==lat[row] && o.longitude_deg==c.nodes[i] || throw(ArgumentError("Observation field or node identity mismatch"))
            end
        end
        new(id,lat,cs,obs)
    end
end
function field_weights(f::EndpointField,p::Real,l::Real)
    p,l=Float64.((p,l));isfinite(p)&&isfinite(l) || throw(ArgumentError("Finite query required"))
    t=(p-f.latitudes[1])/(f.latitudes[2]-f.latitudes[1])
    w=Tuple(a*x for (a,c) in zip((1-t,t),f.charts) for x in cardinal_weights(c,l))
    sum(abs,w)<=AMPLIFICATION_CAP || throw(DomainError((p,l),"Combined weight amplification cap exceeded"))
    w
end
function endpoint(f::EndpointField,p::Real,l::Real)
    w=field_weights(f,p,l)
    (horizontal_m_s=ntuple(k->sum(w[i]*f.observations[i].horizontal_m_s[k] for i in eachindex(w)),2),amplification=sum(abs,w))
end

"""Required physical endpoints and nonzero vertical weights for the query's terrain.
The first lower-table level is computed from terrain, never supplied by a caller.
"""
function required_fields(surface_m::Real,height_m::Real,longitude_deg::Real)
    s,h,l=Float64.((surface_m,height_m,longitude_deg))
    all(isfinite,(s,h,l)) && abs(s)<1e7 && abs(h)<1e7 && 40<=l<=42.5 || throw(ArgumentError("Finite bounded local geometry required"))
    firstlevel=1000floor(Int,s/1000+1.3);clear=h-s
    firstlevel>s+30 && clear>=5 || throw(DomainError(clear,"At least 5 m clearance and valid lower-table geometry required"))
    piece=l<=40.5 ? 1 : 2
    if clear<=30
        a,b,t=FieldID(:clearance,5,piece),FieldID(:clearance,30,piece),(clear-5)/25
    elseif h<=firstlevel
        a,b,t=FieldID(:clearance,30,piece),FieldID(:areoid,firstlevel),(clear-30)/(firstlevel-s-30)
    else
        low=firstlevel+1000floor(Int,(h-firstlevel)/1000)
        a,b,t=FieldID(:areoid,low),FieldID(:areoid,low+1000),(h-low)/1000
    end
    Tuple((field=id,weight=w) for (id,w) in ((a,1-t),(b,t)) if w!=0)
end

"""Immutable local snapshot. Recipe checksum and frozen epoch are required provenance.
Solar offset is hours; the epoch convention is supplied explicitly as text. The
caller must provide a geometry tile with at least 0.25 degrees usable padding.
"""
struct WindSnapshot
    geometry::GeometryTile
    fields::Tuple{Vararg{EndpointField}}
    solar_offset_hours::Float64
    epoch::String
    time_convention::String
    recipe_sha256::String
    function WindSnapshot(geometry::GeometryTile,fields;solar_offset_hours,epoch,time_convention,recipe_sha256)
        fs=Tuple(fields);offset=Float64(solar_offset_hours)
        all(f->f isa EndpointField,fs) && !isempty(fs) && length(unique(f.id for f in fs))==length(fs) || throw(ArgumentError("Unique endpoint fields required"))
        all(x->x>=0.25,geometry.required_padding_deg) || throw(ArgumentError("Quarter-degree slope stencil needs declared geometry padding"))
        isfinite(offset) && 0<=offset<24 && !isempty(epoch) && !isempty(time_convention) && occursin(r"^[0-9a-f]{64}$",recipe_sha256) || throw(ArgumentError("Frozen epoch, time convention, solar offset and recipe checksum required"))
        for level in (5,30),piece in (1,2)
            any(f->f.id==FieldID(:clearance,level,piece),fs) || throw(ArgumentError("Both physical surface endpoints and pieces required"))
        end
        for level in (5,30)
            west=only(f for f in fs if f.id==FieldID(:clearance,level,1))
            east=only(f for f in fs if f.id==FieldID(:clearance,level,2))
            west.latitudes==east.latitudes || throw(ArgumentError("Surface pieces must share fitting latitudes"))
            for row in 1:2
                west.charts[row].hi==east.charts[row].lo==40.5 || throw(ArgumentError("Surface pieces must share the 40.5-degree node"))
                a=west.observations[row*10];b=east.observations[(row-1)*10+1]
                (a.latitude_deg,a.longitude_deg,a.areoid_height_m,a.terrain_height_m,a.horizontal_m_s)==(b.latitude_deg,b.longitude_deg,b.areoid_height_m,b.terrain_height_m,b.horizontal_m_s) || throw(ArgumentError("Shared surface node values disagree"))
            end
        end
        for f in fs,o in f.observations
            native_surface=TerrainGeometry.radial_surface(geometry,o.latitude_deg,o.longitude_deg).terrain_height_m
            abs(native_surface-o.terrain_height_m)<=1e-4 || throw(ArgumentError("Endpoint terrain differs from snapshot geometry"))
        end
        new(geometry,fs,offset,String(epoch),String(time_convention),String(recipe_sha256))
    end
end

# Slope contribution and vertical wind are evaluated using the unclamped
# horizontal wind. The separately supplied sound speed is used afterwards.
function slope_wind(t,p,clearance_km,ew_slope,ns_slope,u,v)
    omega=7.0777e-5;c0=0.07;cosfac=cos(deg2rad(15*(t-15)));depth=3.5+cosfac
    0<clearance_km<=depth || return (0.0,0.0,0.0)
    q=550+150cosfac;xi=clearance_km/depth
    if abs(p)<0.01
        up=(q/c0)*((2+xi*xi)/3-xi)*xi;cross=0.0
    else
        f=2omega*sin(deg2rad(p));wb=2q/(1000depth*abs(f));b=sqrt(1000abs(f)*depth/(2c0))
        den=exp(-4b)-2exp(-2b)*cos(2b)+1
        b1,b2,b3,b4=exp(-b*xi),exp(-b*(2+xi)),exp(-b*(2-xi)),exp(-b*(4-xi))
        up=wb*((b1-b4)*sin(b*xi)+(b2-b3)*sin(b*(2-xi)))/den
        cross=sign(f)*wb*(xi-1+((b1+b4)*cos(b*xi)-(b2+b3)*cos(b*(2-xi)))/den)
    end
    ew=(up*clamp(ew_slope,-.0525,.0525)+cross*clamp(ns_slope,-.0525,.0525))*cosfac
    ns=(up*clamp(ns_slope,-.0525,.0525)-cross*clamp(ew_slope,-.0525,.0525))*sin(deg2rad(15*(t-11)))
    zfac=xi<=.5 ? 1.0 : .5*(1+cos(2pi*(xi-.5)))
    ew,ns,zfac*(ew_slope*(u+ew)+ns_slope*(v+ns))
end

"""Query geodetic latitude, east longitude (degrees), ellipsoid height (metres).
`sound_speed_m_s` is required for native horizontal clamp semantics. Pass
`nothing` explicitly for diagnostic unclamped output. `include_slope=false`
selects background winds. Amplification and provenance remain observable.
"""
function winds(snapshot::WindSnapshot,latitude_deg::Real,longitude_deg::Real,height_m::Real;
               sound_speed_m_s,include_slope::Bool=true)
    c=sound_speed_m_s===nothing ? nothing : Float64(sound_speed_m_s)
    c===nothing || (isfinite(c)&&c>0) || throw(ArgumentError("Positive finite sound speed required"))
    g=TerrainGeometry.geometry(snapshot.geometry,latitude_deg,longitude_deg,height_m)
    p,l=g.planetocentric_latitude_deg,Float64(longitude_deg)
    33.9<=p<=34.6 && 40<=l<=42.5 || throw(DomainError((p,l),"Outside the supported local wind envelope"))
    u=v=0.0;amp=0.0;condition=0.0
    fields=required_fields(g.terrain_height_m,g.areoid_height_m,l)
    for term in fields
        i=findfirst(f->f.id==term.field,snapshot.fields)
        i===nothing && throw(ArgumentError("Missing physical wind field $(term.field)"))
        f=snapshot.fields[i];value=endpoint(f,p,l);u+=term.weight*value.horizontal_m_s[1];v+=term.weight*value.horizontal_m_s[2]
        amp=max(amp,value.amplification);condition=max(condition,maximum(c.condition for c in f.charts))
    end
    background=(u,v,0.0);z=0.0
    if include_slope
        surface(p,l)=TerrainGeometry.radial_surface(snapshot.geometry,p,l).terrain_height_m
        ew=(surface(p,l+.25)-surface(p,l-.25))/(deg2rad(.5)*cos(deg2rad(p))*g.areoid_radius_m)
        ns=(surface(p+.25,l)-surface(p-.25,l))/(deg2rad(.5)*g.areoid_radius_m)
        du,dv,z=slope_wind(mod(snapshot.solar_offset_hours+l/15,24),p,g.radial_clearance_m/1000,ew,ns,u,v)
        u+=du;v+=dv
    end
    unclamped=(u,v,z);limited=c===nothing ? unclamped : (clamp(u,-.7c,.7c),clamp(v,-.7c,.7c),z)
    (wind_m_s=limited,unclamped_m_s=unclamped,background_m_s=background,
     clamp_applied=limited!=unclamped,clamp_enabled=c!==nothing,amplification=amp,condition=condition,
     geometry=g,physical_fields=fields,epoch=snapshot.epoch,time_convention=snapshot.time_convention,recipe_sha256=snapshot.recipe_sha256)
end
end
