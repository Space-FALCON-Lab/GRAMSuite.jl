module ConstructorRecipeTests
using Test, Serialization, StaticArrays, SHA, Libdl

import GRAMSuite
const G = GRAMSuite
const SOURCE = pathof(G)
const EXPECTED_SOURCE_SHA = bytes2hex(sha256(read(SOURCE)))
const LIBRARIES_BEFORE = Set(Libdl.dllist())
# Stand-in for the public Julia binding API. No ccall, dlopen, native library,
# GRAM dataset, SPICE kernel or scientific atmosphere model is used here.
module TraceBinding
const BODY_EARTH = 3
const BODY_MARS = 4
const TRACE = Any[]
const NEXT_ID = Ref(0)
mutable struct Atmosphere
    id::Int
end
record(name, args...; kwargs...) = push!(TRACE, (name=name, args=args, kwargs=(; kwargs...)))
set_library!(path) = (record(:set_library!, path); nothing)
initialize!(path) = (record(:initialize!, path); nothing)
function create_atmosphere(body; kwargs...)
    NEXT_ID[] += 1
    a = Atmosphere(NEXT_ID[])
    record(:create_atmosphere, a, body; kwargs...)
    return a
end
set_start_time!(a; kwargs...) = (record(:set_start_time!, a; kwargs...); nothing)
set_position!(a; kwargs...) = (record(:set_position!, a; kwargs...); nothing)
set_planetary_radii!(a; kwargs...) = (record(:set_planetary_radii!, a; kwargs...); nothing)
set_perturbation_action!(a, flag) = (record(:set_perturbation_action!, a, flag); nothing)
set_seed!(a, seed) = (record(:set_seed!, a, seed); nothing)
set_map_year!(a, year) = (record(:set_map_year!, a, year); nothing)
set_mola_heights!(a, flag) = (record(:set_mola_heights!, a, flag); nothing)
set_perturbation_scales!(a; kwargs...) = (record(:set_perturbation_scales!, a; kwargs...); nothing)
update!(a) = (record(:update!, a); 0)
get_dynamics_state(a) = (density=1.25, temperature=211.0)
get_winds_state(a) = (ewWind=3.0, nsWind=4.0, verticalWind=5.0,
    perturbedEWWind=30.0, perturbedNSWind=40.0, perturbedVerticalWind=50.0)
end

const ROOT = mktempdir(; prefix="gram_standin_root_")
for rel in ("Julia", "Build", "SPICE", "Mars/data", "Earth/data")
    mkpath(joinpath(ROOT, rel))
end
const WRAPPER_MARKER = joinpath(ROOT, "Julia", "GRAM.jl")
write(WRAPPER_MARKER, "# Stand-in marker; never included.\n")
const FAKE_LIBRARY = joinpath(ROOT, "standin-library-marker")
write(FAKE_LIBRARY, "No native code.\n")
# Legacy serialized streams do not carry a library path. Provide only a text
# marker at the resolver's expected path; the stand-in never opens a library.
const AUTO_LIBRARY = G._gram_host_library_path(ROOT)
mkpath(dirname(AUTO_LIBRARY))
write(AUTO_LIBRARY, "No native code; resolver fixture only.\n")
const ORIGINAL_WRAPPER = G._GRAM_WRAPPER[]
const ORIGINAL_WRAPPER_FILE = G._GRAM_WRAPPER_FILE[]
G._GRAM_WRAPPER[] = TraceBinding
G._GRAM_WRAPPER_FILE[] = WRAPPER_MARKER

const EPOCH = G.InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0)
const REFERENCE = (
    initial_time=EPOCH, is_planetocentric=false,
    equatorial_radius_km=3396.190, polar_radius_km=3376.200,
    perturbation_action=false, start_time_scale=1, start_time_frame=0,
    mars_mola_heights=false, mars_map_year=2, seed=4321,
    gram_perturbation_scales=(0.0, 0.0, 0.0, 0.0),
)
function construct(; kwargs...)
    G.GRAMAtmosphereModel(; gram_root_directory=ROOT,
        gram_data_directory=ROOT, spice_directory=joinpath(ROOT, "SPICE"),
        gram_library_path=FAKE_LIBRARY, planet_name="mars", kwargs...)
end
calls(name) = filter(c -> c.name == name, TraceBinding.TRACE)
function lastcall(name)
    cs = calls(name)
    isempty(cs) && error("Missing stand-in call $name")
    return last(cs)
end
function assert_reference_calls(model)
    time = lastcall(:set_start_time!).kwargs
    @test time == (year=2001, month=11, day=6, hour=19, minute=0, seconds=32.0, scale=1, frame=0)
    @test lastcall(:set_planetary_radii!).kwargs == (equatorial=3396.190, polar=3376.200)
    @test lastcall(:set_perturbation_action!).args[2] === false
    @test lastcall(:set_seed!).args[2] == 4321
    @test lastcall(:set_map_year!).args[2] == 2
    @test lastcall(:set_mola_heights!).args[2] === false
    @test lastcall(:set_perturbation_scales!).kwargs ==
        (density_scale=0.0, ew_wind_scale=0.0, ns_wind_scale=0.0, vertical_wind_scale=0.0)
    @test model.is_planetocentric === false
    @test model.offline_surrogate_supported === false
    @test occursin("explicit_native_reference_configuration", model.offline_surrogate_unsupported_reason)
end
function binding_without(name::Symbol)
    mod = Module(gensym(:MissingBindingSetter))
    for n in names(TraceBinding; all=true, imported=false)
        (n == name || n in (:eval, :include, :TraceBinding) || startswith(String(n), "#")) && continue
        isdefined(TraceBinding, n) || continue
        Core.eval(mod, :(const $n = $(getfield(TraceBinding, n))))
    end
    return mod
end
function caught(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

try
results = @testset "native reference API through stand-in binding" begin
    withenv("SPACEAGORA_GRAM_OFFLINE_SURROGATE"=>"off", "SPACEAGORA_GRAM_STATIC_GRID"=>"off",
        "SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_ALL_PLANETS"=>"off", "SPACEAGORA_GRAM_WIND_MODE"=>"nominal") do
        @testset "legacy defaults preserved" begin
            empty!(TraceBinding.TRACE)
            model = construct()
            @test model.is_planetocentric === nothing
            @test model.offline_surrogate_supported
            @test isempty(calls(:set_planetary_radii!))
            @test isempty(calls(:set_perturbation_action!))
            @test lastcall(:set_start_time!).kwargs.scale == 1
            @test lastcall(:set_start_time!).kwargs.frame == 1
            G.point_density_state(model, 150_000.0, deg2rad(45.0), deg2rad(-20.0), 7.0, false)
            @test !haskey(lastcall(:set_position!).kwargs, :is_planetocentric)
        end
        @testset "explicit reference construction and real query paths" begin
            empty!(TraceBinding.TRACE)
            model = construct(; REFERENCE...)
            assert_reference_calls(model)
            expected = (1.25, 211.0, SVector(3.0, 4.0, 5.0))
            @test G.point_density_state(model, 150_000.0, deg2rad(45.0), deg2rad(-20.0), 7.0, false) == expected
            p = lastcall(:set_position!).kwargs
            @test p.height == 150.0
            @test p.latitude ≈ 45.0
            @test p.longitude ≈ -20.0
            @test p.elapsed_time == 7.0
            @test p.is_planetocentric === false
            @test G.density_state(model, 151_000.0, 0.0, 0.0, 9.0, true) == expected
            @test lastcall(:set_position!).kwargs.is_planetocentric === false
            @test lastcall(:set_position!).kwargs.elapsed_time == 9.0
            other = construct(; REFERENCE..., start_time_frame=1, is_planetocentric=true)
            @test model._static_grid_cache_owner !== other._static_grid_cache_owner
            @test lastcall(:set_start_time!).kwargs.frame == 1
            G.point_density_state(other, 150_000.0, 0.0, 0.0, 0.0, false)
            @test lastcall(:set_position!).kwargs.is_planetocentric === true
        end
        @testset "construction recipe copy and serialization" begin
            model = construct(; REFERENCE...)
            empty!(TraceBinding.TRACE)
            copied = deepcopy(model)
            assert_reference_calls(copied)
            @test copied.gram_atmosphere !== model.gram_atmosphere
            @test copied._static_grid_cache_owner !== model._static_grid_cache_owner
            @test copied._constructor_kwargs == model._constructor_kwargs
            @test copied._constructor_kwargs !== model._constructor_kwargs
            G.point_density_state(copied, 150_000.0, 0.0, 0.0, 0.0, false)
            @test lastcall(:set_position!).kwargs.is_planetocentric === false
            io = IOBuffer()
            serialize(io, model)
            seekstart(io)
            empty!(TraceBinding.TRACE)
            restored = deserialize(io)
            assert_reference_calls(restored)
            @test restored.gram_atmosphere !== model.gram_atmosphere
            @test restored._static_grid_cache_owner !== model._static_grid_cache_owner
            @test restored._constructor_kwargs == model._constructor_kwargs
            @test restored._constructor_kwargs !== model._constructor_kwargs
            G.point_density_state(restored, 150_000.0, 0.0, 0.0, 0.0, false)
            @test lastcall(:set_position!).kwargs.is_planetocentric === false
        end
        @testset "legacy positional forms retain unknown recipe" begin
            a = TraceBinding.Atmosphere(-1)
            args = (TraceBinding, a, ROOT, ROOT, joinpath(ROOT, "SPICE"), "mars", EPOCH, true, "")
            models = (G.GRAMAtmosphereModel(args...), G.GRAMAtmosphereModel{Module,typeof(a)}(args...),
                G.GRAMAtmosphereModel(args..., G._GRAMStaticGridCacheOwner()),
                G.GRAMAtmosphereModel{Module,typeof(a)}(args..., G._GRAMStaticGridCacheOwner()))
            for model in models
                @test model._constructor_kwargs === nothing
                @test model.is_planetocentric === nothing
            end
            @test models[1]._static_grid_cache_owner !== models[2]._static_grid_cache_owner
        end
        @testset "invalid reference settings reject before binding calls" begin
            for kwargs in ((equatorial_radius_km=3396.190,), (polar_radius_km=3376.200,),
                (equatorial_radius_km=NaN,polar_radius_km=1.0),
                (equatorial_radius_km=Inf,polar_radius_km=1.0),
                (equatorial_radius_km=0.0,polar_radius_km=0.0),
                (equatorial_radius_km=1.0,polar_radius_km=2.0),
                (start_time_scale=0,), (start_time_frame=2,))
                empty!(TraceBinding.TRACE)
                @test_throws ArgumentError construct(; kwargs...)
                @test isempty(TraceBinding.TRACE)
            end
            @test_throws ArgumentError G.GRAMAtmosphereModel(
                planet_name="earth", equatorial_radius_km=3396.190, polar_radius_km=3376.200)
        end
        @testset "requested missing setters fail clearly before setup" begin
            for (setter, kwargs) in (
                (:set_planetary_radii!, (equatorial_radius_km=3396.190,polar_radius_km=3376.200)),
                (:set_perturbation_action!, (perturbation_action=false,)),
                (:set_position!, (is_planetocentric=false,)),
                (:set_seed!, (is_planetocentric=false,)),
                (:set_map_year!, (mars_map_year=2,)),
                (:set_mola_heights!, (mars_mola_heights=false,)),
                (:set_perturbation_scales!, (gram_perturbation_scales=0.0,)),
            )
                G._GRAM_WRAPPER[] = binding_without(setter)
                empty!(TraceBinding.TRACE)
                err = caught(() -> Base.invokelatest(construct; kwargs...))
                @test err isa ArgumentError
                @test occursin(string(setter), sprint(showerror, err))
                @test isempty(TraceBinding.TRACE)
            end
            G._GRAM_WRAPPER[] = TraceBinding
            # A present function without the required coordinate keyword is
            # also unsupported, rather than silently dropping the request.
            missing_kw = binding_without(:set_position!)
            Core.eval(missing_kw, :(set_position!(a; height, latitude, longitude, elapsed_time) = nothing))
            G._GRAM_WRAPPER[] = missing_kw
            err = caught(() -> Base.invokelatest(construct; is_planetocentric=false))
            @test err isa ArgumentError
            @test occursin("is_planetocentric", sprint(showerror, err))
            bad_action = binding_without(:set_perturbation_action!)
            Core.eval(bad_action, :(set_perturbation_action!(a) = nothing))
            G._GRAM_WRAPPER[] = bad_action
            err = caught(() -> Base.invokelatest(construct; perturbation_action=false))
            @test err isa ArgumentError
            @test occursin("perturbation_action", sprint(showerror, err))
            G._GRAM_WRAPPER[] = TraceBinding
        end
        @testset "serialization schema and old field stream" begin
            old = IOBuffer()
            writer = Serialization.Serializer(old)
            for value in (ROOT, ROOT, joinpath(ROOT, "SPICE"), "mars", EPOCH)
                Serialization.serialize(writer, value)
            end
            seekstart(old)
            restored = Serialization.deserialize(Serialization.Serializer(old), G.GRAMAtmosphereModel)
            @test restored.is_planetocentric === nothing
            @test lastcall(:set_start_time!).kwargs.frame == 1
            bad = IOBuffer()
            Serialization.serialize(Serialization.Serializer(bad), (:GRAMSuite_native_recipe_v99, Dict{Symbol,Any}()))
            seekstart(bad)
            @test_throws ArgumentError Serialization.deserialize(Serialization.Serializer(bad), G.GRAMAtmosphereModel)
        end
    end
    @test Set(Libdl.dllist()) == LIBRARIES_BEFORE
end


finally
    G._GRAM_WRAPPER[] = ORIGINAL_WRAPPER
    G._GRAM_WRAPPER_FILE[] = ORIGINAL_WRAPPER_FILE
    rm(ROOT;recursive=true,force=true)
end

end # module ConstructorRecipeTests
