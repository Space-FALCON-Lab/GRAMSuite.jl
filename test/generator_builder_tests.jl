# These are software tests against a stand-in binding, never native GRAM tests.
using Test, Serialization, TOML, Libdl, Logging

module GeneratorUnderTest
include(joinpath(@__DIR__, "..", "scripts", "build_offline_static_grids.jl"))
end
const GB = GeneratorUnderTest
const GENERATOR_SCRIPT = joinpath(@__DIR__, "..", "scripts", "build_offline_static_grids.jl")
const GENERATOR_LEGACY_SCRIPT = joinpath(@__DIR__, "..", "GRAM Suite 2.0", "simulation", "GRAM", "build_offline_static_grids.jl")

# The stand-in retains a pending/effective MapYear latch so an omitted parameter
# push produces different samples. Its coordinate-dependent fields expose axes.
const GENERATOR_STANDIN = raw"""
module GRAM
export create_atmosphere, set_library!, initialize!, set_start_time!, set_seed!,
    set_perturbation_action!, set_perturbation_scales!, set_planetary_radii!,
    set_map_year!, set_min_max!, set_mola_heights!, set_mgcm_dust_levels!,
    set_dust_storm!, set_dust_density!, set_position!, update!, get_dynamics_state,
    get_winds_state, get_error_message, close!, get_version_string,
    earth_set_merra2_path!
const TRACE = Any[]
const BODY_EARTH=1; const BODY_MARS=2; const BODY_VENUS=3; const BODY_JUPITER=4
const BODY_URANUS=5; const BODY_NEPTUNE=6; const BODY_TITAN=7
mutable struct Atmosphere
    pending::Int
    effective::Int
    position::Any
end
record(name, args...; kwargs...) = push!(TRACE, (name=name, args=args, kwargs=Dict(kwargs)))
function create_atmosphere(body; kwargs...)
    record(:create_atmosphere, body; kwargs...)
    Atmosphere(1, 1, (;height=0.0,latitude=0.0,longitude=0.0,elapsed_time=0.0,is_planetocentric=true))
end
set_library!(path) = record(:set_library!, path)
initialize!(path) = record(:initialize!, path)
for name in (:set_start_time!, :set_seed!, :set_perturbation_action!,
             :set_perturbation_scales!, :set_planetary_radii!, :set_mola_heights!,
             :set_dust_storm!, :set_dust_density!, :earth_set_merra2_path!)
    @eval ($name)(a::Atmosphere, args...; kwargs...) = record($(QuoteNode(name)), args...; kwargs...)
end
function set_map_year!(a, value)
    a.pending=value; record(:set_map_year!, value)
end
function set_min_max!(a, value)
    a.effective=a.pending; record(:set_min_max!, value)
end
function set_mgcm_dust_levels!(a; kwargs...)
    a.effective=a.pending; record(:set_mgcm_dust_levels!; kwargs...)
end
function set_position!(a; kwargs...)
    a.position=(;kwargs...); record(:set_position!; kwargs...)
end
update!(a) = 0
function get_dynamics_state(a)
    p=a.position
    (density=a.effective*1e-8 + p.height*1e-11 + p.latitude*1e-12 + p.longitude*1e-13,
     temperature=180.0+p.height*0.1+p.latitude*0.01)
end
get_winds_state(a) = (ewWind=2.0+a.position.latitude*0.01,
    nsWind=3.0+a.position.longitude*0.01,verticalWind=4.0+a.position.height*0.01)
get_version_string(a) = "stand-in-only\0"
get_error_message() = "stand-in-only"
close!(a) = record(:close!)
end
"""

const GENERATOR_REF_ARGS = ["--planets=mars", "--nalt=2", "--nlat=3", "--nlon=4",
    "--alt-min-km=100", "--alt-max-km=102", "--lat-min-deg=40", "--lat-max-deg=90",
    "--year=2001", "--month=11", "--day=7", "--hour=11", "--minute=51", "--second=4.794789",
    "--elapsed-time-s=0", "--seed=1001", "--is-planetocentric=false",
    "--equatorial-radius-km=3396.190", "--polar-radius-km=3376.200",
    "--start-time-scale=1", "--start-time-frame=0", "--perturbation-action=false",
    "--mars-map-year=2", "--mars-min-max=1", "--mars-mola-heights=false",
    "--mars-dust-mode=native-default"]

function generator_fixture()
    h=[100.0,102.0]; lat=[40.0,90.0]; lon=[0.0,180.0]
    fields=Dict{String,Any}()
    for (key,base,slope) in (("density_kgm3",1e-8,1e-10), ("temperature_K",180.0,0.1),
        ("wind_ew_ms",2.0,0.1),("wind_ns_ms",-3.0,0.2),("wind_up_ms",4.0,0.3))
        fields[key]=[base+slope*(a-100.0) for a in h, b in lat, c in lon]
    end
    Dict{String,Any}("status"=>"ok", "format"=>"spaceagora_gram_static_grid_v1",
        "type"=>"full_grid", "planet"=>"mars", "created_at_utc"=>"2020-01-02T03:04:05",
        "elapsed_s"=>321.0,"elapsed_time_s"=>12.5,"seed"=>7,
        "initial_time"=>Dict("year"=>2001,"month"=>11,"day"=>6,"hour"=>19,"minute"=>0,"second"=>32.0),
        "generation_config"=>Dict("datum"=>"fixture", "nested"=>[1,2]),
        "unknown_metadata"=>Dict("keep"=>[4,5]),
        "grid"=>Dict("alt_km"=>h,"lat_deg"=>lat,"lon_deg"=>lon), "fields"=>fields)
end

function generator_cli(args; root="/nonexistent/generator-test-root", script=GENERATOR_SCRIPT)
    io=IOBuffer()
    command=`$(Base.julia_cmd()) --startup-file=no $script $args`
    process=run(pipeline(ignorestatus(addenv(command,"SPACEAGORA_GRAM_ROOT"=>root)),stdout=io,stderr=io))
    return success(process), String(take!(io))
end

@testset "Offline generator public integration, native-free" begin
    @test !isdefined(GB,:GRAM)
    @test GB.LEGACY_GRAM_ROOT==normpath(joinpath(@__DIR__,"..","GRAM Suite 2.0"))
    # Include the legacy entry in a separate module, with and without an explicit
    # external root. Neither import runs the CLI or leaves an environment change.
    for supplied in (nothing,"/explicit/legacy-test-root")
        withenv("SPACEAGORA_GRAM_ROOT"=>supplied) do
            before=get(ENV,"SPACEAGORA_GRAM_ROOT",nothing)
            legacy=Module(gensym(:LegacyGenerator))
            Base.include(legacy,GENERATOR_LEGACY_SCRIPT)
            @test !isdefined(legacy,:GRAM)
            @test get(ENV,"SPACEAGORA_GRAM_ROOT",nothing)==before
            @test Base.invokelatest(getproperty,legacy,:GRAM_ROOT)==(supplied===nothing ? GB.LEGACY_GRAM_ROOT : supplied)
        end
    end
    @test !any(occursin("libgram",lowercase(path)) for path in Libdl.dllist())
    @testset "Strict CLI and preflight" begin
        @test GB.prepare_build(GENERATOR_REF_ARGS).cfg.start_time_frame == 0
        @test GB.prepare_build(GENERATOR_REF_ARGS).native_required
        @test GB.prepare_build(["--build-grid=false","--build-surrogate=true","--surrogate-source=grid"]).native_required == false
        @test GB.parse_cli(["--strict"])["strict"] == "true"
        @test GB.parse_cli(["--strict", "false"])["strict"] == "false"
        @test GB.parse_cli(["--second", "4.794789"])["second"] == "4.794789"
        for args in (["--typo=1"],["--nlat=3","--nlat=4"],["--nlat"],["--nlat="],
            ["nlat=3"],["--strict=maybe"],["--build-grid=maybe"],["--surrogate-only=maybe"],
            ["--is-planetocentric=maybe"],["--perturbation-action=maybe"])
            @test_throws ArgumentError GB.parse_cli(args)
        end
        for args in (["--nalt=1"],["--nlon=1"],["--nlat=1"],["--lat-min-deg=91"],
            ["--alt-min-km=NaN"],["--alt-max-km=Inf"],["--elapsed-time-s=Inf"],
            ["--year=2001","--month=2","--day=29"],["--hour=24"],["--second=60"],
            ["--seed=-1"],["--seed=2147483648"],["--start-time-frame=2"],
            ["--start-time-scale=0"],["--mars-map-year=3"],["--mars-min-max=0"],
            ["--equatorial-radius-km=3396.19"],["--surrogate-lon-step-deg=360"],
            ["--surrogate-lat-step-deg=Inf"],["--surrogate-adaptive-pilot-alt-step-km=Inf"],
            ["--planets=venus","--equatorial-radius-km=3","--polar-radius-km=2"],
            ["--planets=venus","--mars-map-year=2"],
            ["--build-surrogate=true","--surrogate-alt-bands=0:Inf:1"])
            @test_throws Exception GB.prepare_build(args)
        end
        for value in ("40", "40,40", "90,40", "40,NaN,90", "40,91", "40,,90")
            @test_throws ArgumentError GB.prepare_build(["--lat-nodes-deg=$value"])
        end
        for key in ("nlat=4","lat-min-deg=40","lat-max-deg=90","surrogate-lat-step-deg=2.5",
            "surrogate-adaptive-pilot-lat-step-deg=5")
            @test_throws ArgumentError GB.prepare_build(["--lat-nodes-deg=40,65,90","--$key"])
        end
        cfg=GB.prepare_build(["--lat-nodes-deg=40,42.8225,60.2795,77.6361,90"]).cfg
        @test cfg.nlat==5
        @test GB.latitude_axis(cfg;count=3)==[40,42.8225,60.2795,77.6361,90]
        mktempdir() do root
            explicit=joinpath(root,"explicit")
            fallback=joinpath(root,"Earth","data","MERRA2data")
            mkpath(joinpath(fallback,"All Mean"))
            write(joinpath(fallback,"All Mean","MERRA2All_11.bin"),"fixture-only")
            @test_throws ArgumentError GB.resolve_earth_merra2_path(root,11;override=explicit)
            withenv("EARTH_MERRA2_PATH"=>nothing) do
                @test GB.resolve_earth_merra2_path(root,11)==fallback
                @test GB.resolve_earth_merra2_path(root,12)===nothing
            end
            mkpath(joinpath(explicit,"All Mean"))
            write(joinpath(explicit,"All Mean","MERRA2All_11.bin"),"fixture-only")
            @test GB.resolve_earth_merra2_path(root,11;override=explicit)==explicit
            @test_throws ArgumentError GB.resolve_earth_merra2_path(root,12;override=explicit)
            @test_throws ArgumentError GB.resolve_earth_merra2_path(root,0;override=explicit)
        end
        @test !isdefined(GB,:GRAM)
        good,log=generator_cli(["--help"])
        @test good
        @test occursin("SPACEAGORA_GRAM_ROOT",log)
        @test occursin("--mars-min-max",log)
        good,log=generator_cli(["--help"];script=GENERATOR_LEGACY_SCRIPT)
        @test good
        @test occursin("Advanced frozen-atmosphere grid generator",log)
        @test occursin("--mars-min-max",log)
        bad,log=generator_cli(["--nalt=2","--nalt=3"])
        @test !bad
        @test occursin("Duplicate option",log)
        @test !occursin("Native GRAM binding not found",log)
        bad,log=generator_cli(["--lat-min-deg=91"])
        @test !bad
        @test occursin("Latitude bounds",log)
        @test !occursin("Native GRAM binding not found",log)
        bad,log=generator_cli(GENERATOR_REF_ARGS)
        @test !bad
        @test occursin("Set SPACEAGORA_GRAM_ROOT",log)
    end
    @testset "Resampling preserves source and rejects unsupported bounds" begin
        source=generator_fixture(); before=deepcopy(source)
        cfg=GB.prepare_build(GENERATOR_REF_ARGS).cfg
        target=(precision=Float64,alt_nodes=[100.0,101.0,102.0],alt_bands=[(100.0,102.0,1.0)],
            lat_step_deg=25.0,lon_step_deg=90.0,alt_mode="bands")
        product=with_logger(NullLogger()) do
            GB.build_surrogate_from_grid_payload(source,"mars",cfg,target)
        end
        @test size(product["fields"]["density_kgm3"])==(3,3,4)
        @test product["fields"]["density_kgm3"][:,1,1]≈[1e-8,1.01e-8,1.02e-8]
        @test product["elapsed_time_s"]==12.5
        @test product["seed"]==7
        @test product["initial_time"]==source["initial_time"]
        @test product["generation_config"]==source["generation_config"]
        @test product["unknown_metadata"]==source["unknown_metadata"]
        @test product["build_metadata"]["requested_config"]["elapsed_time_s"]==0
        @test product["build_metadata"]["source_metadata"]["created_at_utc"]==source["created_at_utc"]
        @test source==before
        product["unknown_metadata"]["keep"][1]=99
        product["generation_config"]["nested"][1]=99
        @test source==before
        for cfgbad in (merge(cfg,(lat_min_deg=39.0,)),merge(cfg,(lat_max_deg=90.1,)))
            @test_throws ArgumentError GB.build_surrogate_from_grid_payload(source,"mars",cfgbad,target)
        end
        for targetbad in (merge(target,(alt_nodes=[99.0,102.0],)),merge(target,(alt_nodes=[100.0,103.0],)),
            merge(target,(alt_nodes=[102.0,100.0],)),merge(target,(alt_nodes=[100.0,NaN],)))
            @test_throws ArgumentError GB.build_surrogate_from_grid_payload(source,"mars",cfg,targetbad)
        end
        @test_throws ArgumentError GB.build_surrogate_from_grid_payload(source,"venus",cfg,target)
        for key in ("density_kgm3","temperature_K","wind_ew_ms","wind_ns_ms","wind_up_ms")
            broken=deepcopy(source);broken["fields"][key][1]=NaN
            @test_throws ArgumentError GB.build_surrogate_from_grid_payload(broken,"mars",cfg,target)
        end
        @test GB.lookup_payload_state(source,101.0,65.0,-1.0)==GB.lookup_payload_state(source,101.0,65.0,359.0)
        mktempdir() do dir
            serialize(joinpath(dir,"mars_grid.jls"),source)
            args=["--planets=mars","--out-dir=$dir","--build-grid=false","--build-surrogate=true",
                "--surrogate-source=grid","--strict=true","--alt-min-km=100","--alt-max-km=102",
                "--lat-min-deg=40","--lat-max-deg=90","--surrogate-lat-step-deg=25",
                "--surrogate-lon-step-deg=90","--surrogate-alt-bands=100:102:1"]
            good,log=generator_cli(args)
            @test good
            replay=deserialize(joinpath(dir,"surrogates","mars_surrogate.jls"))
            @test replay["fields"]["density_kgm3"]≈Float32.(product["fields"]["density_kgm3"])
            @test replay["seed"]==7
            @test TOML.parsefile(joinpath(dir,"grid_build_summary.toml"))["planet"][1]["surrogate_source"]=="grid"
            good,log=generator_cli(args;script=GENERATOR_LEGACY_SCRIPT)
            @test good
            legacy_replay=deserialize(joinpath(dir,"surrogates","mars_surrogate.jls"))
            @test legacy_replay["fields"]==replay["fields"]
            @test legacy_replay["grid"]==replay["grid"]
            @test legacy_replay["generation_config"]==replay["generation_config"]
            @test legacy_replay["seed"]==replay["seed"]
        end
        # Adaptive profiles inherit the entire source axis. Reject an altitude
        # subset instead of emitting full-range data under narrower request metadata.
        mktempdir() do dir
            wide=generator_fixture()
            wide["grid"]["alt_km"]=[0.0,300.0]
            serialize(joinpath(dir,"mars_grid.jls"),wide)
            common=["--planets=mars","--out-dir=$dir","--build-grid=false",
                "--build-surrogate=true","--surrogate-source=grid","--strict=true",
                "--surrogate-adaptive-alt=true","--surrogate-adaptive-min-step-km=50",
                "--surrogate-adaptive-max-step-km=50","--lat-nodes-deg=40,65,90",
                "--surrogate-lon-step-deg=90"]
            good,log=generator_cli(vcat(common,["--alt-min-km=100","--alt-max-km=200"]))
            @test !good
            @test occursin("Adaptive altitude nodes span [0.0, 300.0] km instead of requested [100.0, 200.0] km",log)
            @test !occursin("Native GRAM binding not found",log)
            rejected=deserialize(joinpath(dir,"surrogates","mars_surrogate.jls"))
            @test rejected["status"]=="error"
            @test !haskey(rejected,"fields")
            @test deserialize(joinpath(dir,"mars_grid.jls"))==wide
            good,log=generator_cli(vcat(common,["--alt-min-km=0","--alt-max-km=300"]))
            @test good
            accepted=deserialize(joinpath(dir,"surrogates","mars_surrogate.jls"))
            @test accepted["grid"]["alt_km"]==collect(0.0:50.0:300.0)
            @test accepted["grid"]["lat_deg"]==[40.0,65.0,90.0]
            @test accepted["surrogate"]["alt_mode"]=="adaptive_profile"
            @test accepted["build_metadata"]["requested_config"]["alt_min_km"]==0.0
            @test accepted["build_metadata"]["requested_config"]["alt_max_km"]==300.0
            @test accepted["seed"]==wide["seed"]
        end
        mktempdir() do dir
            good,log=generator_cli(["--planets=mars","--out-dir=$dir","--build-grid=false",
                "--build-surrogate=true","--surrogate-source=grid"])
            @test !good
            @test occursin("Grid source unavailable",log)
            @test !occursin("Native GRAM binding not found",log)
            @test !occursin("falling back",log)
        end
        @test !isdefined(GB,:GRAM)
    end
    @testset "Actual source using a stand-in binding" begin
        Core.eval(GB,Meta.parse(GENERATOR_STANDIN))
        Core.eval(GB,:(using .GRAM))
        cfg=GB.prepare_build(GENERATOR_REF_ARGS).cfg
        trace=GB.GRAM.TRACE
        calls(name)=filter(x->x.name==name,trace)
        target=(precision=Float64,alt_nodes=[100.0,102.0],alt_bands=[(100.0,102.0,2.0)],
            lat_step_deg=25.0,lon_step_deg=90.0,alt_mode="bands")
        for route in (:full,:direct,:pilot)
            empty!(trace)
            output = if route==:full
                Base.invokelatest(GB.build_planet_grid,"unused","mars",cfg)
            elseif route==:direct
                Base.invokelatest(GB.build_surrogate_direct,"unused","mars",cfg,target)
            else
                Base.invokelatest(GB.density_profile_direct,"unused","mars",cfg;
                    pilot_alt_step_km=2.0,pilot_lat_step_deg=25.0,pilot_lon_step_deg=90.0)
            end
            names=getproperty.(trace,:name)
            @test findfirst(==(:set_map_year!),names)<findfirst(==(:set_min_max!),names)<findfirst(==(:set_position!),names)
            @test calls(:set_map_year!)[1].args==(2,)
            @test calls(:set_min_max!)[1].args==(1,)
            @test calls(:set_start_time!)[1].kwargs[:frame]==0
            @test calls(:set_start_time!)[1].kwargs[:scale]==1
            @test calls(:set_start_time!)[1].kwargs[:seconds]==4.794789
            @test calls(:set_planetary_radii!)[1].kwargs==Dict(:equatorial=>3396.190,:polar=>3376.200)
            @test calls(:set_mola_heights!)[1].args==(false,)
            @test calls(:set_perturbation_action!)[1].args==(false,)
            @test all(==(0.0),values(calls(:set_perturbation_scales!)[1].kwargs))
            @test isempty(calls(:set_mgcm_dust_levels!))
            @test isempty(calls(:set_dust_storm!))
            @test isempty(calls(:set_dust_density!))
            @test length(calls(:set_position!))==24
            @test all(p->p.kwargs[:is_planetocentric]===false,calls(:set_position!))
            @test all(p->p.kwargs[:elapsed_time]==0.0,calls(:set_position!))
            @test sort(unique(p.kwargs[:latitude] for p in calls(:set_position!)))==[40.0,65.0,90.0]
            @test length(calls(:close!))==1
            if route != :pilot
                @test output["generation_config"]["native_model_version"]=="stand-in-only"
                @test output["generation_config"]["native_model_version_status"]=="native_getter"
                @test output["generation_config"]["mars_map_year"]==2
                @test output["generation_config"]["mars_min_max"]==1
                @test output["fields"]["density_kgm3"][1,1,1]≈2e-8+100e-11+40e-12
            end
        end
        for nodes in ([40.0,42.8225,60.2795,77.6361,90.0],[40.0,65.0,90.0])
            explicit=merge(cfg,(lat_nodes_deg=nodes,))
            product=Base.invokelatest(GB.build_planet_grid,"unused","mars",explicit)
            @test product["grid"]["lat_deg"]==nodes
            @test product["latitude_sampling"]["lat_nodes_count"]==length(nodes)
            @test product["grid"]["lat_deg"]!==nodes
        end
        empty!(trace)
        @test_throws ArgumentError Base.invokelatest(GB.build_planet_grid,"unused","mars",merge(cfg,(lat_nodes_deg=[90.0,40.0],)))
        @test isempty(trace)
        # Omit the parameter push in the modelled latch as a negative control.
        pending=Base.invokelatest(GB.GRAM.create_atmosphere,2)
        Base.invokelatest(GB.GRAM.set_map_year!,pending,2)
        @test pending.pending==2 && pending.effective==1
        Base.invokelatest(GB.configure_planet_specific!,"mars",pending,cfg)
        @test pending.effective==2
        # Exercise lazy binding import and the actual script main in a fresh Julia.
        mktempdir() do root
            mkpath(joinpath(root,"Julia"));write(joinpath(root,"Julia","GRAM.jl"),GENERATOR_STANDIN)
            out=joinpath(root,"output")
            good,log=generator_cli(vcat(GENERATOR_REF_ARGS,["--out-dir=$out","--strict=true"]);root=root)
            @test good
            if good
                product=deserialize(joinpath(out,"mars_grid.jls"))
                @test size(product["fields"]["density_kgm3"])==(2,3,4)
                @test product["generation_config"]["is_planetocentric"]===false
                @test product["generation_config"]["mars_min_max"]==1
                summary=TOML.parsefile(joinpath(out,"grid_build_summary.toml"))
                @test summary["planet"][1]["native_model_version"]=="stand-in-only"
                @test summary["planet"][1]["native_model_version_status"]=="native_getter"
            else
                @info "Stand-in CLI failure" log
            end
        end
    end
    @test !any(occursin("libgram",lowercase(path)) for path in Libdl.dllist())
end
