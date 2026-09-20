module TerrainGeometryTests
using Test
using GRAMSuite.TerrainGeometry

const E=Ellipsoid(3396190.,3376200.)
const D=ServedDomain((34.,35.),(40.,42.5),(-5000.,41250.))
const lat=[33.25,33.75,34.25,34.75,35.25]
const lon=[39.488+0.5j for j in 0:7]
const a=fill(3.39e6,5,8)
const t=zeros(5,8)
const corners=[radial_position(E,p,h).planetocentric_latitude_deg
               for p in D.geodetic_latitude_deg for h in D.ellipsoid_height_m]

@testset "Reviewer boundary cases carried into port" begin
    for side in (:west,:east,:south,:north),n in (128,129)
        la,lo=copy(lat),copy(lon)
        if side==:west
            lo=collect(range(40.0-n*eps(43.),42.988,length=8))
        elseif side==:east
            lo=collect(range(39.488,42.5+n*eps(43.),length=8))
        elseif side==:south
            la=collect(range(minimum(corners)-n*eps(35.25),35.25,length=5))
        else
            la=collect(range(33.25,maximum(corners)+n*eps(35.25),length=5))
        end
        request=side in (:west,:east) ? (.25,64eps(43.)) : (64eps(35.25),.25)
        if n==128
            @test_throws ArgumentError GeometryTile(la,lo,a,t,E,D;minimum_padding_deg=request)
        else
            tile=GeometryTile(la,lo,a,t,E,D;minimum_padding_deg=request)
            @test tile isa GeometryTile
            for p in D.geodetic_latitude_deg,l in D.east_longitude_deg,h in D.ellipsoid_height_m
                @test isfinite(geometry(tile,p,l,h).radial_clearance_m)
            end
        end
    end
    d10=ServedDomain((10.,11.),(40.,42.5),(-5000.,41250.))
    tile=GeometryTile([9.25,9.75,10.25,10.75,11.25],lon,a,t,E,d10;minimum_padding_deg=(.25,.25))
    c=padding_certificate(tile)
    @test c.numerical_floor_deg==(64eps(11.25),64eps(42.988))
    @test c.numerical_floor_deg[1] != c.numerical_floor_deg[2]
    for i in 1:4
        floor=c.numerical_floor_deg[i<=2 ? 1 : 2]
        raw=c.actual_south_north_west_east_deg[i]
        guard=c.guarded_south_north_west_east_deg[i]
        @test guard==prevfloat(raw-floor)
        @test BigFloat(guard)<=BigFloat(raw)-BigFloat(floor)
    end
    for n in (0,64,65)
        lo=copy(lon);lo[1]=39.75-n*eps(43.);lo[end]=42.75+n*eps(43.)
        if n<65
            @test_throws ArgumentError GeometryTile(lat,lo,a,t,E,D;minimum_padding_deg=(.25,.25))
        else
            @test GeometryTile(lat,lo,a,t,E,D;minimum_padding_deg=(.25,.25)) isa GeometryTile
        end
    end
    @test ServedDomain((-89.,-88.),(40.,42.5),(-100000.,0.)) isa ServedDomain
    @test ServedDomain((88.,89.),(40.,42.5),(-100000.,0.)) isa ServedDomain
    @test_throws ArgumentError ServedDomain((prevfloat(-89.),-88.),(40.,42.5),(-100000.,0.))
    @test_throws ArgumentError ServedDomain((88.,nextfloat(89.)),(40.,42.5),(-100000.,0.))
    @test_throws ArgumentError ServedDomain((-1.,1.),(40.,42.5),(prevfloat(-100000.),0.))
    @test ServedDomain((-1.,1.),(0.,359.9),(-5000.,0.)) isa ServedDomain
    @test_throws ArgumentError ServedDomain((-1.,1.),(0.,360.),(-5000.,0.))
    for bounds in ((-89.,-88.),(88.,89.))
        domain=ServedDomain(bounds,(40.,42.5),(-100000.,100000.))
        cs=[radial_position(E,p,h).planetocentric_latitude_deg for p in bounds for h in domain.ellipsoid_height_m]
        tile=GeometryTile([minimum(cs)-.3,maximum(cs)+.3],[39.5,43.],fill(3.39e6,2,2),zeros(2,2),E,domain;minimum_padding_deg=(.25,.25))
        for p in bounds,l in domain.east_longitude_deg,h in domain.ellipsoid_height_m
            @test isfinite(geometry(tile,p,l,h).radial_clearance_m)
        end
    end
end

deeply_immutable(x) = !ismutabletype(typeof(x)) &&
    all(deeply_immutable(getfield(x,i)) for i in 1:fieldcount(typeof(x)))

@testset "Immutable storage and native-free geometry" begin
    # Cross terms exercise latitude-major indexing and all four weights.
    radius(p,l)=3.39e6+2p+3l+0.1p*l
    terrain(p,l)=-100+p-l+0.03p*l
    la,lo=copy(lat),copy(lon)
    ar=[radius(p,l) for p in la,l in lo]
    tr=[terrain(p,l) for p in la,l in lo]
    tile=GeometryTile(la,lo,ar,tr,E,D;minimum_padding_deg=(.25,.25))
    before=geometry(tile,34.5,41.,0.)
    certificate=padding_certificate(tile)
    @test deeply_immutable(tile)
    @test size(tile.areoid_radius_m)==(5,8)
    @test tile.areoid_radius_m[2,3]==ar[2,3]
    @test collect(tile.latitude_deg)==la
    @test Matrix(tile.terrain_height_m)==tr

    # Caller-owned arrays can change after construction without changing the tile.
    la[1]=NaN;lo[end]=NaN;ar[2,3]=NaN;tr[1,1]=NaN
    @test geometry(tile,34.5,41.,0.)==before
    @test padding_certificate(tile)==certificate
    for storage in (tile.latitude_deg,tile.longitude_deg,tile.areoid_radius_m,tile.terrain_height_m)
        snapshot=copy(storage)
        @test_throws Exception setindex!(storage,NaN,1)
        @test_throws Exception fill!(storage,NaN)
        @test_throws Exception setindex!(view(storage,:),NaN,1)
        @test_throws Exception setindex!(getfield(storage,:data),NaN,1)
        @test copy(storage)==snapshot
        snapshot[1]=NaN
        @test all(isfinite,storage)
    end
    @test_throws Exception setproperty!(tile,:latitude_deg,copy(lat))
    @test geometry(tile,34.5,41.,0.)==before
    @test padding_certificate(tile)==certificate

    for p in range(34.,35.,length=7),l in range(40.,42.5,length=9),h in (-5000.,0.,41250.)
        q=geometry(tile,p,l,h)
        @test isapprox(q.areoid_radius_m,radius(q.planetocentric_latitude_deg,l);atol=1e-8,rtol=0)
        @test isapprox(q.terrain_height_m,terrain(q.planetocentric_latitude_deg,l);atol=1e-10,rtol=0)
        @test q.radial_clearance_m==q.radius_m-q.areoid_radius_m-q.terrain_height_m
    end
    queries=[(p,l,h) for p in (34.,34.5,35.) for l in (40.,41.,42.5) for h in (-5000.,0.,41250.)]
    expected=[geometry(tile,q...) for q in queries]
    tasks=[Threads.@spawn([geometry(tile,q...) for q in queries]) for _ in 1:8]
    @test all(fetch(task)==expected for task in tasks)

    for query in ((prevfloat(34.),41.,0.),(nextfloat(35.),41.,0.),
                  (34.5,prevfloat(40.),0.),(34.5,nextfloat(42.5),0.),
                  (34.5,41.,prevfloat(-5000.)),(34.5,41.,nextfloat(41250.)),
                  (NaN,41.,0.),(34.5,Inf,0.),(34.5,41.,NaN))
        @test_throws DomainError geometry(tile,query...)
    end
    @test_throws ArgumentError GeometryTile(reverse(lat),lon,a,t,E,D;minimum_padding_deg=(.25,.25))
    @test_throws DimensionMismatch GeometryTile(lat,lon,a[1:2,:],t,E,D;minimum_padding_deg=(.25,.25))
    @test_throws UndefKeywordError GeometryTile(lat,lon,a,t,E,D)
    for padding in ((0.,.25),(.25,0.),(-1.,.25),(NaN,.25),(.25,Inf),(.25,),(.25,.25,.25))
        @test_throws ArgumentError GeometryTile(lat,lon,a,t,E,D;minimum_padding_deg=padding)
    end
    @test_throws ArgumentError GeometryTile(lat,lon,fill(NaN,5,8),t,E,D;minimum_padding_deg=(.25,.25))
    @test_throws ArgumentError GeometryTile(lat,lon,fill(-1.,5,8),t,E,D;minimum_padding_deg=(.25,.25))
    @test_throws ArgumentError GeometryTile(lat,lon,a,fill(Inf,5,8),E,D;minimum_padding_deg=(.25,.25))
    @test_throws ArgumentError GeometryTile(lat,lon,a,t,E,D;minimum_padding_deg=(.5,.5))
    @test_throws DomainError radial_position(E,34.,-E.polar_m^2/E.equatorial_m)
end
end # module
