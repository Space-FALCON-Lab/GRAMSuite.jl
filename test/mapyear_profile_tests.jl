module MapYearProfileTests
using Test, Serialization, StaticArrays, SHA, Libdl

import GRAMSuite
const G = GRAMSuite
const SOURCE = pathof(G)
const EXPECTED_SOURCE_SHA = bytes2hex(sha256(read(SOURCE)))
const LIBRARIES_BEFORE = Set(Libdl.dllist())
# This stand-in models only the pending/effective configuration distinction
# reported by the retained native probes. It does not model atmospheric physics.
# State markers returned as density/temperature make a missing or premature push
# observable through the real wrapper query path, not merely through a call log.
module LatchedBinding
const BODY_MARS = 4
const TRACE = Any[]
const NEXT_ID = Ref(0)
mutable struct Atmosphere
    id::Int
    pending_map_year::Int
    effective_map_year::Int
    min_max::Union{Nothing,Int}
end
record(name, a::Atmosphere, args...; kwargs...) =
    push!(TRACE, (name=name, id=a.id, args=args, kwargs=(; kwargs...)))
set_library!(path) = nothing
initialize!(path) = nothing
function create_atmosphere(body; kwargs...)
    NEXT_ID[] += 1
    a = Atmosphere(NEXT_ID[], 1, 1, nothing)
    record(:create_atmosphere, a, body; kwargs...)
    return a
end
set_start_time!(a; kwargs...) = (record(:set_start_time!, a; kwargs...); nothing)
set_position!(a; kwargs...) = (record(:set_position!, a; kwargs...); nothing)
set_seed!(a, value) = (record(:set_seed!, a, value); nothing)
set_perturbation_action!(a, value) = (record(:set_perturbation_action!, a, value); nothing)
set_perturbation_scales!(a; kwargs...) = (record(:set_perturbation_scales!, a; kwargs...); nothing)
set_planetary_radii!(a; kwargs...) = (record(:set_planetary_radii!, a; kwargs...); nothing)
set_mola_heights!(a, value) = (record(:set_mola_heights!, a, value); nothing)
function set_map_year!(a, year::Integer)
    record(:set_map_year!, a, year)
    a.pending_map_year = Int(year)
    return nothing
end
function set_min_max!(a, value::Integer)
    record(:set_min_max!, a, value)
    a.min_max = Int(value)
    a.effective_map_year = a.pending_map_year
    return nothing
end
update!(a) = (record(:update!, a); 0)
get_dynamics_state(a) = (density=Float64(a.effective_map_year), temperature=100.0+a.effective_map_year)
get_winds_state(a) = (ewWind=0.0, nsWind=0.0, verticalWind=0.0,
    perturbedEWWind=0.0, perturbedNSWind=0.0, perturbedVerticalWind=0.0)
end

const ROOT = mktempdir(; prefix="gram_standin_root_")
for rel in ("Julia", "Build", "SPICE", "Mars/data")
    mkpath(joinpath(ROOT, rel))
end
const MARKER = joinpath(ROOT, "Julia", "GRAM.jl")
const LIBRARY_MARKER = joinpath(ROOT, "not-a-native-library")
write(MARKER, "# Never included: stand-in binding installed explicitly.\n")
write(LIBRARY_MARKER, "No native code.\n")
const ORIGINAL_WRAPPER = G._GRAM_WRAPPER[]
const ORIGINAL_WRAPPER_FILE = G._GRAM_WRAPPER_FILE[]
G._GRAM_WRAPPER[] = LatchedBinding
G._GRAM_WRAPPER_FILE[] = MARKER

# Fixture expresses the named Odyssey profile. This is not a new wrapper
# default or an attempt to make all generic map-year-only requests safe.
const PREVIOUS_PROFILE = (
    initial_time=G.InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0),
    is_planetocentric=false, equatorial_radius_km=3396.190, polar_radius_km=3376.200,
    perturbation_action=false, start_time_scale=1, start_time_frame=0,
    mars_mola_heights=false, mars_map_year=2, seed=1001,
    gram_perturbation_scales=(0.0,0.0,0.0,0.0),
)
const PROFILE = merge(PREVIOUS_PROFILE, (mars_min_max=1,))
construct(options) = G.GRAMAtmosphereModel(;
    gram_root_directory=ROOT, gram_data_directory=ROOT,
    spice_directory=joinpath(ROOT, "SPICE"), gram_library_path=LIBRARY_MARKER,
    planet_name="mars", options...)
query(model) = G.point_density_state(model, 80_000.0, deg2rad(-25.0), deg2rad(240.0), 0.0, false)
trace(model) = filter(t -> t.id == model.gram_atmosphere.id, LatchedBinding.TRACE)
function assert_effective_profile(model)
    a = model.gram_atmosphere
    @test a.pending_map_year == 2
    @test a.effective_map_year == 2
    @test a.min_max == 1
    calls = trace(model)
    maps = findall(t -> t.name == :set_map_year!, calls)
    pushes = findall(t -> t.name == :set_min_max!, calls)
    @test length(maps) == 1
    @test length(pushes) == 1
    @test only(maps) < only(pushes)
    @test calls[only(maps)].args == (2,)
    @test calls[only(pushes)].args == (1,)
    @test model._constructor_kwargs[:mars_min_max] == 1
    @test model._constructor_kwargs[:mars_map_year] == 2
    @test query(model) == (2.0,102.0,SVector(0.0,0.0,0.0))
    @test G.density_state(model, 80_000.0, 0.0, 0.0, 30.0, true)[1] == 2.0
    queries = findall(t -> t.name == :set_position!, trace(model))
    @test !isempty(queries) && all(i -> i > only(pushes), queries)
end

try
results = @testset "explicit Odyssey MapYear application profile" begin
    withenv("SPACEAGORA_GRAM_STATIC_GRID"=>"off", "SPACEAGORA_GRAM_OFFLINE_SURROGATE"=>"off",
        "SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_ALL_PLANETS"=>"off", "SPACEAGORA_GRAM_WIND_MODE"=>"nominal") do
        @testset "negative control retains pending map-year-only state" begin
            previous = construct(PREVIOUS_PROFILE)
            @test previous.gram_atmosphere.pending_map_year == 2
            @test previous.gram_atmosphere.effective_map_year == 1
            @test previous.gram_atmosphere.min_max === nothing
            @test isempty(filter(t -> t.name == :set_min_max!, trace(previous)))
            @test previous._constructor_kwargs[:mars_min_max] === nothing
            for _ in 1:3
                @test query(previous)[1] == 1.0
            end
            @test G.density_state(previous, 80_000.0, 0.0, 0.0, 0.0, true)[1] == 1.0
            @test previous.gram_atmosphere.effective_map_year == 1
            generic = construct((mars_map_year=2,))
            @test query(generic)[1] == 1.0
            @test generic.gram_atmosphere.min_max === nothing
        end
        @testset "wrong setter order is discriminating" begin
            a = LatchedBinding.create_atmosphere(LatchedBinding.BODY_MARS)
            LatchedBinding.set_min_max!(a,1)
            LatchedBinding.set_map_year!(a,2)
            LatchedBinding.update!(a)
            @test a.pending_map_year == 2
            @test a.effective_map_year == 1
            @test LatchedBinding.get_dynamics_state(a).density == 1.0
            LatchedBinding.set_min_max!(a,1)
            @test a.effective_map_year == 2
            @test LatchedBinding.get_dynamics_state(a).density == 2.0
        end
        @testset "explicit profile, copy and serialization preserve application" begin
            model = construct(PROFILE)
            assert_effective_profile(model)
            copied = deepcopy(model)
            assert_effective_profile(copied)
            @test copied.gram_atmosphere !== model.gram_atmosphere
            @test copied._constructor_kwargs == model._constructor_kwargs
            @test copied._constructor_kwargs !== model._constructor_kwargs
            @test copied._static_grid_cache_owner !== model._static_grid_cache_owner
            io = IOBuffer(); serialize(io,model); seekstart(io)
            restored = deserialize(io)
            assert_effective_profile(restored)
            @test restored.gram_atmosphere !== model.gram_atmosphere
            @test restored._constructor_kwargs == model._constructor_kwargs
            @test restored._constructor_kwargs !== model._constructor_kwargs
            @test restored._static_grid_cache_owner !== model._static_grid_cache_owner
        end
        @testset "missing explicitly requested setter fails" begin
            missing = Module(:MissingMinMaxBinding)
            for name in names(LatchedBinding; all=true, imported=false)
                (name in (:eval,:include,:LatchedBinding,:set_min_max!) || startswith(String(name),"#")) && continue
                isdefined(LatchedBinding,name) || continue
                Core.eval(missing, :(const $name = $(getfield(LatchedBinding,name))))
            end
            G._GRAM_WRAPPER[] = missing
            before = LatchedBinding.NEXT_ID[]
            err = try
                Base.invokelatest(construct,PROFILE)
                nothing
            catch exception
                exception
            end
            @test err isa ArgumentError
            @test occursin("mars_min_max",sprint(showerror,err))
            @test occursin("set_min_max!",sprint(showerror,err))
            @test LatchedBinding.NEXT_ID[] == before
            G._GRAM_WRAPPER[] = LatchedBinding
        end
    end
    @test bytes2hex(sha256(read(SOURCE))) == EXPECTED_SOURCE_SHA
    @test Set(Libdl.dllist()) == LIBRARIES_BEFORE
end


finally
    G._GRAM_WRAPPER[] = ORIGINAL_WRAPPER
    G._GRAM_WRAPPER_FILE[] = ORIGINAL_WRAPPER_FILE
    rm(ROOT;recursive=true,force=true)
end

end # module MapYearProfileTests
