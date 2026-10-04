using LinearAlgebra
const MW=GRAMSuite.MarsWindEndpoints

function synthetic_wind_snapshot()
    G=GRAMSuite.TerrainGeometry
    lat=[33.25,33.75,34.25,34.75,35.25];lon=collect(39.488:.5:42.988)
    terrain=[-2000+4*(p-34.2)+7*(l-41) for p in lat,l in lon]
    tile=G.GeometryTile(lat,lon,fill(3390000.,length(lat),length(lon)),terrain,G.Ellipsoid(3396190,3376200),G.ServedDomain((34.19,34.93),(40,42.5),(-5000,10000));minimum_padding_deg=(.25,.25))
    lats=(34.15,34.35);fields=MW.EndpointField[]
    for id in (MW.FieldID(:areoid,-1000),MW.FieldID(:areoid,0),MW.FieldID(:areoid,1000),MW.FieldID(:areoid,2000),MW.FieldID(:areoid,3000),[MW.FieldID(:clearance,z,p) for z in (5,30) for p in (1,2)]...)
        lo,hi=id.datum===:areoid ? (40.,42.5) : id.piece==1 ? (40.,40.5) : (40.5,42.5)
        charts=(MW.HarmonicChart(lo,hi,id.datum===:areoid ? 5 : 10),MW.HarmonicChart(lo,hi,id.datum===:areoid ? 5 : 10))
        observations=MW.Observation[]
        for (p,c) in zip(lats,charts),l in c.nodes
            s=G.radial_surface(tile,p,l).terrain_height_m;h=id.datum===:areoid ? id.level_m : s+id.level_m
            # A known field in the declared harmonic space with linear latitude.
            value=(10+2sin(deg2rad(l))+.1p+id.level_m*.001, -4+cos(2deg2rad(l))-.2p)
            push!(observations,MW.Observation(id,p,l,h,s,value))
        end
        push!(fields,MW.EndpointField(id,lats,charts,observations))
    end
    snapshot=MW.WindSnapshot(tile,fields;solar_offset_hours=12.48,epoch="synthetic fixed epoch",time_convention="synthetic fixed solar offset",recipe_sha256=repeat("a",64))
    snapshot,fields,terrain
end

@testset "Immutable physical endpoint winds, synthetic data only" begin
    s,caller_fields,caller_terrain=synthetic_wind_snapshot()
    west=only(f for f in s.fields if f.id==MW.FieldID(:clearance,5,1))
    obs=collect(west.observations);o=obs[10]
    obs[10]=MW.Observation(o.field,o.latitude_deg,o.longitude_deg,o.areoid_height_m,o.terrain_height_m,(o.horizontal_m_s[1]+1,o.horizontal_m_s[2]))
    wrong=MW.EndpointField(west.id,west.latitudes,west.charts,obs)
    badfields=[f.id==west.id ? wrong : f for f in s.fields]
    @test_throws ArgumentError MW.WindSnapshot(s.geometry,badfields;solar_offset_hours=12.,epoch="synthetic",time_convention="UTC",recipe_sha256=repeat("a",64))
    @test_throws ArgumentError MW.WindSnapshot(s.geometry,[s.fields...,first(s.fields)];solar_offset_hours=12.,epoch="synthetic",time_convention="UTC",recipe_sha256=repeat("a",64))
    for f in s.fields
        for c in f.charts
            @test isfinite(c.condition) && c.condition<1000
            for (i,node) in enumerate(c.nodes)
                w=MW.cardinal_weights(c,node)
                @test maximum(abs.(collect(w).-Float64.(1:c.size .== i)))<1e-11
            end
        end
        p,l=34.58,40.35
        expected=(10+2sin(deg2rad(l))+.1p+f.id.level_m*.001,-4+cos(2deg2rad(l))-.2p)
        if f.id.datum===:areoid || f.id.piece==1
            @test maximum(abs.(MW.endpoint(f,p,l).horizontal_m_s.-expected))<1e-9
        end
    end
    # An independent high-precision cardinal solve in the ill-conditioned raw
    # trig basis verifies stable Float64 chart weights during extrapolation.
    setprecision(BigFloat,256) do
        c=MW.HarmonicChart(40.15,41.42,5);l=41.43
        plain(x)=BigFloat[1,sin(deg2rad(x)),cos(deg2rad(x)),sin(2deg2rad(x)),cos(2deg2rad(x))]
        a=reduce(vcat,permutedims.(plain.(BigFloat.(c.nodes))))
        w=transpose(a)\plain(BigFloat(l))
        @test maximum(abs.(w.-collect(MW.cardinal_weights(c,l))))<big"1e-10"
    end
    for z in (5,30)
        terms=MW.required_fields(-2000,-2000+z,40.5)
        @test length(terms)==1 && terms[1].field==MW.FieldID(:clearance,z,1)
    end
    @test only(MW.required_fields(-2000,-1000,40.5)).field==MW.FieldID(:areoid,-1000)
    @test only(MW.required_fields(-2000,0,40.5)).field==MW.FieldID(:areoid,0)
    @test_throws DomainError MW.required_fields(-2000,-1995.1,41)
    @test_throws ArgumentError MW.required_fields(NaN,0,41)
    @test_throws ArgumentError MW.FieldID(:ordinal,1000)
    @test_throws ArgumentError MW.FieldID(:areoid,1)
    @test_throws ArgumentError MW.HarmonicChart(40,50,5)
    @test_throws ArgumentError MW.HarmonicChart(40,41,6)
    @test_throws DomainError MW.cardinal_weights(MW.HarmonicChart(40,41,10),42)
    f=first(s.fields)
    @test_throws DomainError MW.field_weights(f,100.,40.4)
    @test_throws ArgumentError MW.EndpointField(f.id,f.latitudes,f.charts,f.observations[1:end-1])
    @test_throws ArgumentError MW.EndpointField(MW.FieldID(:areoid,0),f.latitudes,f.charts,f.observations)
    @test_throws ArgumentError MW.Observation(f.id,34.2,41,-1000,-500,(1,2))
    @test_throws ArgumentError MW.Observation(f.id,34.2,41,-1000,-2000,(NaN,2))
    @test_throws MethodError setindex!(f.observations,first(f.observations),1)
    @test_throws ErrorException setfield!(first(f.observations),:terrain_height_m,0.)
    @test_throws ErrorException setfield!(f.charts[1],:condition,1.)
    @test_throws MethodError setindex!(f.charts[1].inverse[1],0.,1)
    # Find a height with 17.3 m clearance using the accepted geometry path.
    low,high=-5000.,10000.
    for _ in 1:60
        mid=(low+high)/2;g=GRAMSuite.TerrainGeometry.geometry(s.geometry,34.5,41.1,mid)
        if g.radial_clearance_m<17.3;low=mid;else;high=mid;end
    end
    h=(low+high)/2
    original=MW.winds(s,34.5,41.1,h;sound_speed_m_s=nothing)
    limited=MW.winds(s,34.5,41.1,h;sound_speed_m_s=1.)
    @test limited.clamp_applied
    @test maximum(abs,limited.wind_m_s[1:2])<=.7
    @test limited.wind_m_s[3]==original.wind_m_s[3] != 0
    @test limited.unclamped_m_s==original.wind_m_s
    @test original.amplification<=10 && original.condition<1000
    @test original.epoch==s.epoch && original.recipe_sha256==s.recipe_sha256
    empty!(caller_fields);fill!(caller_terrain,0)
    @test MW.winds(s,34.5,41.1,h;sound_speed_m_s=nothing)==original
    @test all(fetch(t)==original for t in [Threads.@spawn(MW.winds(s,34.5,41.1,h;sound_speed_m_s=nothing)) for _ in 1:16])
    @test_throws ArgumentError MW.winds(s,34.5,41.1,h;sound_speed_m_s=0.)
    @test_throws DomainError MW.winds(s,34.,41.1,h;sound_speed_m_s=nothing)
    @test_throws DomainError MW.winds(s,34.5,43.,h;sound_speed_m_s=nothing)
    @test_throws ArgumentError MW.WindSnapshot(s.geometry,s.fields;solar_offset_hours=12.,epoch="",time_convention="UTC",recipe_sha256=repeat("a",64))
    reduced=MW.WindSnapshot(s.geometry,filter(f->f.id.datum===:clearance,s.fields);solar_offset_hours=12.,epoch="synthetic",time_convention="UTC",recipe_sha256=repeat("a",64))
    @test_throws ArgumentError MW.winds(reduced,34.5,41.1,h+500;sound_speed_m_s=nothing)
end
