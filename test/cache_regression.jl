module CacheRegressionTests
using Test
using Serialization
using SHA
using Libdl

import GRAMSuite
const G = GRAMSuite
const CACHE_SOURCE = pathof(G)
const LIBRARIES_BEFORE = Set(Libdl.dllist())

mutable struct FakeAtmosphere
    marker::Float64
    seed::Int
end
FakeAtmosphere(marker) = FakeAtmosphere(Float64(marker), 1001)

const BUILD_COUNTS = IdDict{Any,Int}()
const BUILT_MODELS = Any[]

# Test only the production cache lifecycle. This more-specific builder cannot
# call native GRAM and creates tiny grids with independently identifiable data.
function G._gram_static_grid_build(
    model::G.GRAMAtmosphereModel{M,FakeAtmosphere}, key::G.GRAMStaticGridKey;
    lock_obj=nothing
) where {M}
    owner = getfield(model, :_static_grid_cache_owner)
    BUILD_COUNTS[owner] = get(BUILD_COUNTS, owner, 0) + 1
    push!(BUILT_MODELS, model)
    dims = (key.n_alt, key.n_lat, key.n_lon)
    marker = model.initial_time.year + model.gram_atmosphere.marker
    effective_wind = key.include_wind && G._gram_wind_mode() == :perturbed ? 23.0 : 7.0
    return G.GRAMStaticGrid(
        key,
        collect(range(key.alt_min_m, key.alt_max_m; length=key.n_alt)),
        collect(range(-pi/2, pi/2; length=key.n_lat)),
        collect(range(0.0, 2pi; length=key.n_lon + 1))[1:end-1],
        fill(marker, dims), fill(200.0, dims),
        fill(effective_wind, dims), zeros(dims), zeros(dims),
    )
end

function make_model(; year=2000, marker=1.0, planet="mars", handle=FakeAtmosphere(marker), explicit=false)
    arguments = (nothing, handle, "", "", "", planet,
                 G.InitialTime(year=Int32(year)), true, "")
    return explicit ? G.GRAMAtmosphereModel{Nothing,FakeAtmosphere}(arguments...) :
                      G.GRAMAtmosphereModel(arguments...)
end

build_count(model) = get(BUILD_COUNTS, getfield(model, :_static_grid_cache_owner), 0)
cache_entry_count() = sum(owner -> length(owner.entries), keys(G._GRAM_STATIC_GRID_CACHE); init=0)

@noinline function cache_owners_then_release(count)
    owners = WeakRef[]
    arrays = WeakRef[]
    for index in 1:count
        model = make_model(marker=index)
        grid = G._gram_static_grid_get_or_build!(model, true)
        push!(owners, WeakRef(getfield(model, :_static_grid_cache_owner)))
        for array in (grid.rho, grid.T, grid.wind_e, grid.wind_n, grid.wind_u)
            push!(arrays, WeakRef(array))
        end
    end
    # Test instrumentation must not retain models or tokens during this proof.
    empty!(BUILT_MODELS)
    empty!(BUILD_COUNTS)
    return owners, arrays
end

@noinline function keep_model_alive_through_gc()
    model = make_model(marker=75)
    grid = G._gram_static_grid_get_or_build!(model, true)
    owner_ref = WeakRef(getfield(model, :_static_grid_cache_owner))
    array_ref = WeakRef(grid.rho)
    grid = nothing
    empty!(BUILT_MODELS)
    empty!(BUILD_COUNTS)
    GC.gc(true)
    @test owner_ref.value !== nothing
    @test array_ref.value !== nothing
    # GRAMStaticGrid is immutable and may be reboxed; its mutable array identity
    # is the stable observation. Instrumentation was reset before collection.
    @test G._gram_static_grid_get_or_build!(model, true).rho === array_ref.value
    @test build_count(model) == 0
    empty!(BUILT_MODELS)
    empty!(BUILD_COUNTS)
    return owner_ref, array_ref
end

function argument_error_message(f)
    err = try
        f()
        nothing
    catch error
        error
    end
    @test err isa ArgumentError
    return err isa Exception ? sprint(showerror, err) : ""
end

# The actual keyword constructor, deepcopy and custom deserialization are tested
# with a deliberately substituted Julia backend. Its methods never load a
# library or read GRAM assets; all directories and placeholder files are local
# test fixtures. This is not native construction validation.
module FakeNativeBackend
const BODY_EARTH = 1
const BODY_MARS = 2
const CREATED = Ref(0)
set_library!(path) = nothing
initialize!(path) = nothing
function create_atmosphere(body; data_path=nothing)
    CREATED[] += 1
    return getfield(parentmodule(@__MODULE__), :FakeAtmosphere)(7.0)
end
set_start_time!(handle; kwargs...) = nothing
set_seed!(handle, seed) = (handle.seed = seed; nothing)
end

@testset "Per-model runtime static-grid ownership" begin
    @test G._GRAM_WRAPPER[] === nothing
    withenv(
        "GRAM_ROOT" => nothing, "GRAM_LIB" => nothing,
        "SPACEAGORA_GRAM_STATIC_GRID" => "off",
        "SPACEAGORA_GRAM_STATIC_GRID_ALT_MIN_M" => "0",
        "SPACEAGORA_GRAM_STATIC_GRID_ALT_MAX_M" => "1000",
        "SPACEAGORA_GRAM_STATIC_GRID_NALT" => "3",
        "SPACEAGORA_GRAM_STATIC_GRID_NLAT" => "3",
        "SPACEAGORA_GRAM_STATIC_GRID_NLON" => "8",
        "SPACEAGORA_GRAM_STATIC_GRID_ELAPSED_TIME_S" => "0",
        "SPACEAGORA_GRAM_STATIC_GRID_WIND" => "true",
        "SPACEAGORA_GRAM_WIND_MODE" => "nominal",
        "SPACEAGORA_GRAM_STATIC_GRID_PLANETS" => nothing,
        "SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_ALL_PLANETS" => "off",
        "SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_STRICT" => "true",
    ) do
        G.clear_gram_static_grid_cache!()
        @testset "Epoch separation, same-model hits and constructor compatibility" begin
            first_model = make_model()
            later_model = make_model(year=2010)
            @test G._gram_static_grid_key(first_model, true) == G._gram_static_grid_key(later_model, true)
            first_grid = G._gram_static_grid_get_or_build!(first_model, true)
            later_grid = G._gram_static_grid_get_or_build!(later_model, true)
            @test first_grid !== later_grid
            @test first_grid.rho[1] == 2001.0
            @test later_grid.rho[1] == 2011.0
            @test G._gram_static_grid_get_or_build!(first_model, true) === first_grid
            @test G._gram_static_grid_get_or_build!(later_model, true) === later_grid
            @test build_count(first_model) == build_count(later_model) == 1

            shared_handle = FakeAtmosphere(2.0)
            inferred = make_model(handle=shared_handle)
            parameterized = make_model(handle=shared_handle, explicit=true)
            @test typeof(inferred) == typeof(parameterized)
            @test inferred.gram_atmosphere === parameterized.gram_atmosphere
            @test getfield(inferred, :_static_grid_cache_owner) !== getfield(parameterized, :_static_grid_cache_owner)
            @test G._gram_static_grid_get_or_build!(inferred, true) !== G._gram_static_grid_get_or_build!(parameterized, true)
            @test build_count(inferred) == build_count(parameterized) == 1
        end

        @testset "Same epoch, incompatible settings, mutation and selective invalidation" begin
            G.clear_gram_static_grid_cache!()
            first_model = make_model(marker=10.0)
            other_model = make_model(marker=100.0)
            first_grid = G._gram_static_grid_get_or_build!(first_model, true)
            other_grid = G._gram_static_grid_get_or_build!(other_model, true)
            @test first_grid !== other_grid
            @test first_grid.rho[1] == 2010.0
            @test other_grid.rho[1] == 2100.0
            withenv("SPACEAGORA_GRAM_STATIC_GRID_NALT" => "4") do
                additional = G._gram_static_grid_get_or_build!(first_model, true)
                @test additional !== first_grid
                @test size(additional.rho) == (4, 3, 8)
            end
            @test cache_entry_count() == 3
            first_model.gram_atmosphere.marker = 50.0
            # Raw native state changes are deliberately not auto-detected.
            @test G._gram_static_grid_get_or_build!(first_model, true) === first_grid
            @test G.clear_gram_static_grid_cache!(first_model) === nothing
            @test cache_entry_count() == 1
            @test G._gram_static_grid_get_or_build!(other_model, true) === other_grid
            @test build_count(other_model) == 1
            replacement = G._gram_static_grid_get_or_build!(first_model, true)
            @test replacement !== first_grid
            @test replacement.rho[1] == 2050.0
            @test build_count(first_model) == 3
            @test G.clear_gram_static_grid_cache!() === nothing
            @test isempty(G._GRAM_STATIC_GRID_CACHE)
            @test G._gram_static_grid_get_or_build!(first_model, true) !== replacement
            @test G._gram_static_grid_get_or_build!(other_model, true) !== other_grid
            @test build_count(first_model) == 4 && build_count(other_model) == 2
            @test cache_entry_count() == 2
            G.clear_gram_static_grid_cache!()
        end

        @testset "Effective wind-mode keys distinguish sampled policies" begin
            G.clear_gram_static_grid_cache!()
            model = make_model()
            nominal = G._gram_static_grid_get_or_build!(model, true)
            @test nominal.wind_e[1] == 7.0
            withenv("SPACEAGORA_GRAM_WIND_MODE" => "perturbed") do
                perturbed = G._gram_static_grid_get_or_build!(model, true)
                @test perturbed !== nominal
                @test perturbed.wind_e[1] == 23.0
                @test build_count(model) == 2
                withenv("SPACEAGORA_GRAM_WIND_MODE" => "auto") do
                    @test G._gram_static_grid_get_or_build!(model, true) === perturbed
                end
            end
            withenv("SPACEAGORA_GRAM_WIND_MODE" => "mean") do
                @test G._gram_static_grid_get_or_build!(model, true) === nominal
            end
            @test build_count(model) == 2
            withenv("SPACEAGORA_GRAM_STATIC_GRID_WIND" => "false") do
                wind_false_nominal = G._gram_static_grid_get_or_build!(model, true)
                @test wind_false_nominal.wind_e[1] == 7.0
                withenv("SPACEAGORA_GRAM_WIND_MODE" => "perturbed") do
                    @test G._gram_static_grid_get_or_build!(model, true) === wind_false_nominal
                end
            end
            @test build_count(model) == 3
            withenv("SPACEAGORA_GRAM_WIND_MODE" => "invalid") do
                @test_throws ArgumentError G._gram_static_grid_get_or_build!(model, true)
            end
            @test build_count(model) == 3
            @test cache_entry_count() == 3
            G.clear_gram_static_grid_cache!(model)
            @test cache_entry_count() == 0
        end

        @testset "Weak ownership reclaims dropped grids while preserving live models" begin
            G.clear_gram_static_grid_cache!()
            empty!(BUILT_MODELS)
            empty!(BUILD_COUNTS)
            live_owner, live_array = keep_model_alive_through_gc()
            GC.gc(true)
            GC.gc(true)
            @test live_owner.value === nothing
            @test live_array.value === nothing
            @test isempty(G._GRAM_STATIC_GRID_CACHE)
            owners, arrays = cache_owners_then_release(50)
            @test length(owners) == 50
            @test length(arrays) == 250
            GC.gc(true)
            GC.gc(true)
            @test all(ref -> ref.value === nothing, owners)
            @test all(ref -> ref.value === nothing, arrays)
            @test isempty(G._GRAM_STATIC_GRID_CACHE)
            println("gc_lifecycle: released owners=", count(ref -> ref.value === nothing, owners),
                    "/", length(owners), ", released grid arrays=", count(ref -> ref.value === nothing, arrays),
                    "/", length(arrays))
        end

        @testset "Precompute warms supplied instances and rejects temporary-model promises" begin
            G.clear_gram_static_grid_cache!()
            empty!(BUILT_MODELS)
            mars = make_model(planet="mars")
            earth = make_model(planet="earth")
            @test G.precompute_gram_static_grids!(mars) === nothing
            @test length(BUILT_MODELS) == 1 && BUILT_MODELS[1] === mars
            @test build_count(mars) == 1
            @test G.precompute_gram_static_grids!(mars; planets=["mars"]) === nothing
            @test build_count(mars) == 1
            before = length(BUILT_MODELS)
            message = argument_error_message(() -> G.precompute_gram_static_grids!(mars; planets=["mars", "earth"]))
            @test occursin("Single-model precomputation", message)
            @test occursin("precompute_gram_static_grids!([", message)
            @test length(BUILT_MODELS) == before
            withenv("SPACEAGORA_GRAM_STATIC_GRID_PLANETS" => "earth") do
                @test_throws ArgumentError G.precompute_gram_static_grids!(mars)
                @test length(BUILT_MODELS) == before
            end
            @test G.precompute_gram_static_grids!([mars, earth]) === nothing
            @test build_count(mars) == build_count(earth) == 1
            @test length(BUILT_MODELS) == 2 && BUILT_MODELS[2] === earth
            @test G._gram_static_grid_get_or_build!(mars, true) === G._gram_static_grid_get_or_build!(mars, true)
            @test G._gram_static_grid_get_or_build!(earth, true) === G._gram_static_grid_get_or_build!(earth, true)
            @test_throws ArgumentError G.precompute_gram_static_grids!(G.GRAMAtmosphereModel[])
            mixed_case = make_model(planet="MaRs")
            @test G.precompute_gram_static_grids!(mixed_case; planets=["mars"]) === nothing
            @test build_count(mixed_case) == 1
            withenv("SPACEAGORA_GRAM_STATIC_GRID" => "on",
                    "SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_ALL_PLANETS" => "on") do
                message = argument_error_message(() -> G.GRAMAtmosphereModel(; gram_root_directory="does-not-exist"))
                @test occursin("SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_ALL_PLANETS", message)
                @test occursin("precompute_gram_static_grids!([", message)
                @test G._GRAM_WRAPPER[] === nothing
            end
        end

        @testset "Keyword reconstruction, deepcopy and deserialization get fresh owners" begin
            mktempdir() do root
                for folder in ("Build/lib", "Julia", "SPICE", "Earth/data", "Mars/data")
                    mkpath(joinpath(root, folder))
                end
                wrapper_file = joinpath(root, "Julia", "GRAM.jl")
                write(wrapper_file, "# Test placeholder; never included.\n")
                library = joinpath(root, "Build", "lib", G._gram_host_library_filename())
                write(library, "Test placeholder; never loaded as a library.\n")
                G._GRAM_WRAPPER[] = FakeNativeBackend
                G._GRAM_WRAPPER_FILE[] = wrapper_file
                try
                    model = G.GRAMAtmosphereModel(; gram_root_directory=root,
                        gram_data_directory=root, spice_directory=joinpath(root, "SPICE"),
                        planet_name="mars", initial_time=G.InitialTime(year=2015))
                    original_grid = G._gram_static_grid_get_or_build!(model, true)
                    copied = deepcopy(model)
                    @test copied.initial_time == model.initial_time
                    @test copied.gram_atmosphere !== model.gram_atmosphere
                    @test getfield(copied, :_static_grid_cache_owner) !== getfield(model, :_static_grid_cache_owner)
                    @test G._gram_static_grid_get_or_build!(copied, true) !== original_grid

                    io = IOBuffer()
                    serialize(io, model)
                    seekstart(io)
                    restored = deserialize(io)
                    @test restored.initial_time == model.initial_time
                    @test restored.gram_atmosphere !== model.gram_atmosphere
                    @test getfield(restored, :_static_grid_cache_owner) !== getfield(model, :_static_grid_cache_owner)
                    @test G._gram_static_grid_get_or_build!(restored, true) !== original_grid
                    @test FakeNativeBackend.CREATED[] == 3
                finally
                    G._GRAM_WRAPPER[] = nothing
                    G._GRAM_WRAPPER_FILE[] = ""
                end
            end
        end
        G.clear_gram_static_grid_cache!()
        @test isempty(G._GRAM_STATIC_GRID_CACHE)
    end
    @test G._GRAM_WRAPPER[] === nothing
    @test isempty(filter(path -> occursin(r"(?i)(libgram|libcspice)", path),
                         collect(setdiff(Set(Libdl.dllist()), LIBRARIES_BEFORE))))
end

println("source_sha256=", bytes2hex(sha256(read(CACHE_SOURCE))))
println("native_backend_used=false; constructor reconstruction tested with explicit Julia test double")

end # module CacheRegressionTests
