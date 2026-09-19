using Test
using SHA
using TOML
using Serialization

include(joinpath(@__DIR__, "..", "scripts", "grid_generation_recipe.jl"))
const GRecipe = GridGenerationRecipe
const RECIPE = joinpath(@__DIR__, "..", "recipes", "odyssey_p20_frozen_v1.toml")

module RecipeBuilderContract
include(joinpath(@__DIR__, "..", "scripts", "build_offline_static_grids.jl"))
end

function assert_manifest(directory)
    lines = readlines(joinpath(directory, "SHA256SUMS"))
    names = String[]
    for line in lines
        sha, name = split(line, "  "; limit=2)
        @test sha == bytes2hex(open(sha256, joinpath(directory, name)))
        @test name != "SHA256SUMS"
        push!(names, name)
    end
    actual = [replace(relpath(joinpath(d, f), directory), '\\' => '/') for (d, _, fs) in walkdir(directory) for f in fs if joinpath(d, f) != joinpath(directory, "SHA256SUMS")]
    @test sort(names) == sort(actual)
    @test length(names) == length(unique(names))
end

@testset "Recipe is the retained Odyssey fine configuration" begin
    r = GRecipe.read_recipe(RECIPE)
    opts = GRecipe.builder_options(r, "outputs")
    @test opts["planets"] == "mars"
    @test opts["nalt"] == "161"
    @test opts["nlat"] == "51"
    @test opts["nlon"] == "144"
    @test opts["year"] == "2001"
    @test opts["month"] == "11"
    @test opts["day"] == "7"
    @test opts["hour"] == "11"
    @test opts["minute"] == "51"
    @test opts["second"] == "4.794789"
    @test opts["start-time-frame"] == "0"
    @test opts["start-time-scale"] == "1"
    @test opts["elapsed-time-s"] == "0.0"
    @test opts["is-planetocentric"] == "false"
    @test opts["mars-mola-heights"] == "false"
    @test opts["mars-map-year"] == "2"
    @test opts["mars-min-max"] == "1"
    @test opts["mars-dust-mode"] == "native-default"
    @test opts["equatorial-radius-km"] == "3396.19"
    @test opts["polar-radius-km"] == "3376.2"
    @test opts["perturbation-action"] == "false"
    @test opts["seed"] == "1001"
    @test opts["strict"] == "true"
    @test opts["build-grid"] == "true"
    @test opts["build-surrogate"] == "false"
    @test r["grid"]["field_type"] == "Float64"
    lock = GRecipe.read_input_lock(r, RECIPE)
    @test length(lock.lock["planet_data"]) == 9
    @test length(lock.lock["spice_directory_inventory"]) == 30
    @test lock.lock["planet"] == "mars"
    @test sort(GRecipe.argument_vector(opts)) == GRecipe.argument_vector(opts)
    plan = RecipeBuilderContract.prepare_build(GRecipe.argument_vector(opts))
    @test plan.native_required
    @test plan.build_grid && !plan.build_surrogate
    @test (plan.cfg.nalt, plan.cfg.nlat, plan.cfg.nlon) == (161, 51, 144)
    @test plan.cfg.start_time_frame == 0
    @test plan.cfg.mars_map_year == 2 && plan.cfg.mars_min_max == 1
    @test !isdefined(RecipeBuilderContract, :GRAM)
end

@testset "Invalid scientific settings fail before execution" begin
    base = GRecipe.read_recipe(RECIPE)
    mutations = [
        r -> (r["unknown"] = true),
        r -> delete!(r["time"], "epoch_utc"),
        r -> (r["schema_version"] = true),
        r -> (r["schema_version"] = 1.0),
        r -> (r["time"]["scale_code"] = true),
        r -> (r["time"]["epoch_utc"] = "2001-02-30T00:00:00Z"),
        r -> (r["time"]["epoch_utc"] = "2001-11-07T11:51:60Z"),
        r -> (r["time"]["epoch_utc"] = "2001-11-07T11:51:04+00:00"),
        r -> (r["time"]["epoch_utc"] = "2001-11-07T11:51:04.1234567890Z"),
        r -> (r["time"]["elapsed_seconds"] = 1.0),
        r -> (r["time"]["frame_code"] = 1),
        r -> (r["time"]["scale"] = "ET"),
        r -> (r["coordinates"]["latitude"] = "planetocentric"),
        r -> (r["coordinates"]["height"] = "MOLA"),
        r -> (r["coordinates"]["radii_mode"] = "native-default"),
        r -> delete!(r["coordinates"], "polar_radius_km"),
        r -> (r["coordinates"]["polar_radius_km"] = Inf),
        r -> (r["coordinates"]["polar_radius_km"] = 4000.0),
        r -> (r["atmosphere"]["seed"] = true),
        r -> (r["atmosphere"]["perturbation_action"] = true),
        r -> (r["atmosphere"]["perturbation_scales"] = [0.0, 0.1, 0.0, 0.0]),
        r -> (r["atmosphere"]["wind_product"] = "perturbed"),
        r -> delete!(r["mars"], "min_max"),
        r -> (r["mars"]["min_max"] = 0),
        r -> (r["mars"]["map_year"] = 3),
        r -> (r["mars"]["mola_heights"] = true),
        r -> (r["grid"]["altitude_step_km"] = 3.0),
        r -> (r["grid"]["altitude_step_km"] = NaN),
        r -> (r["grid"]["latitude_bounds_deg"] = [40.0, 91.0]),
        r -> (r["grid"]["longitude_bounds_deg"] = [10.0, 360.0]),
        r -> (r["grid"]["longitude_step_deg"] = 360.0),
        r -> (r["grid"]["field_type"] = "Float32"),
        r -> (r["provenance"]["native_library_sha256"] = "not-a-hash"),
        r -> (r["provenance"]["input_lock_file"] = "../outside.toml"),
        r -> delete!(r["provenance"], "input_lock_sha256"),
        r -> (r["provenance"]["unprobed_defaults"] = String[])]
    for mutate in mutations
        recipe = deepcopy(base); mutate(recipe)
        @test_throws ArgumentError GRecipe.validate_recipe(recipe)
    end
    for planet in filter(!=("mars"), collect(GRecipe.PLANETS))
        r = deepcopy(base)
        r["planet"] = planet; delete!(r, "mars")
        r["coordinates"]["radii_mode"] = "native-default"
        delete!(r["coordinates"], "equatorial_radius_km")
        delete!(r["coordinates"], "polar_radius_km")
        @test GRecipe.validate_recipe(r) === r # Syntax support is not a supported scientific preset.
        @test !haskey(GRecipe.builder_options(r, "outputs"), "mars-map-year")
    end
    @test_throws ArgumentError GRecipe.main(["--recipe=a", "--recipe=b"])
    @test_throws ArgumentError GRecipe.main(["--unexpected=a"])
end

# A real Julia subprocess emits a tiny synthetic builder payload. It never loads
# the native binding, library, native model or supplied data files.
const STANDIN = raw"""
using Serialization, TOML
opts = Dict(split(a[3:end], '='; limit=2) for a in ARGS)
@assert ENV["GRAM_LIB"] == opts["lib"]
@assert realpath(pwd()) == realpath(dirname(opts["out-dir"]))
@assert realpath(Base.active_project()) == realpath(joinpath(pwd(), "child_project", "Project.toml"))
@assert LOAD_PATH == ["@", "@stdlib"]
if isfile(joinpath(ENV["SPACEAGORA_GRAM_ROOT"], "FAIL")); error("deliberate stand-in failure"); end
if isfile(joinpath(ENV["SPACEAGORA_GRAM_ROOT"], "MUTATE")); write(opts["lib"], "modified during stand-in execution"); end
out = opts["out-dir"]; mkpath(out)
nalt, nlat, nlon = (parse(Int, opts[k]) for k in ("nalt", "nlat", "nlon"))
alt = collect(range(parse(Float64, opts["alt-min-km"]), parse(Float64, opts["alt-max-km"]); length=nalt))
lat = collect(range(parse(Float64, opts["lat-min-deg"]), parse(Float64, opts["lat-max-deg"]); length=nlat))
lon = collect(range(0.0, 360.0; length=nlon+1))[1:end-1]
applied = Dict{String, Any}("start_time_scale" => 1, "start_time_frame" => 0, "elapsed_time_s" => 0.0, "seed" => parse(Int, opts["seed"]), "is_planetocentric" => false, "perturbation_action" => false, "perturbation_scales" => zeros(4), "initial_time" => Dict{String, Any}(k => (k == "second" ? parse(Float64, opts[k]) : parse(Int, opts[k])) for k in ("year", "month", "day", "hour", "minute", "second")), "native_model_version" => "synthetic stand-in; not GRAM")
for key in ("map-year", "min-max"); applied["mars_" * replace(key, "-" => "_")] = parse(Int, opts["mars-" * key]); end
applied["mars_mola_heights"] = false; applied["mars_dust_mode"] = opts["mars-dust-mode"]
for key in ("equatorial-radius-km", "polar-radius-km"); applied[replace(key, "-" => "_")] = parse(Float64, opts[key]); end
p = Dict("status" => "ok", "planet" => "mars", "grid" => Dict("alt_km" => alt, "lat_deg" => lat, "lon_deg" => lon), "fields" => Dict(k => fill(1.0, nalt, nlat, nlon) for k in ("density_kgm3", "temperature_K", "wind_ew_ms", "wind_ns_ms", "wind_up_ms")), "generation_config" => applied)
if isfile(joinpath(ENV["SPACEAGORA_GRAM_ROOT"], "WRONG")); p["generation_config"]["mars_min_max"] = 0; end
open(joinpath(out, "mars_grid.jls"), "w") do io; serialize(io, p); end
open(opts["summary-file"], "w") do io; TOML.print(io, Dict("stand_in" => true)); end
println("synthetic subprocess completed without loading GRAM")
"""

@testset "Atomic runs, input inventories and real stand-in subprocess" begin
    mktempdir() do tmp
        pkg = joinpath(tmp, "source package")
        mkpath(joinpath(pkg, "scripts"))
        for name in ("generate_grid.jl", "grid_generation_recipe.jl")
            cp(joinpath(@__DIR__, "..", "scripts", name), joinpath(pkg, "scripts", name))
        end
        write(joinpath(pkg, "scripts", "build_offline_static_grids.jl"), STANDIN)
        native = joinpath(tmp, "fake native")
        for dir in ("Julia", "Mars/data", "SPICE/lsk", "Earth/data")
            mkpath(joinpath(native, dir))
        end
        write(joinpath(native, "Julia", "GRAM.jl"), "error(\"this must never be loaded\")")
        write(joinpath(native, "Mars", "data", "table.bin"), "selected Mars input")
        write(joinpath(native, "Earth", "data", "table.bin"), "unselected Earth input")
        write(joinpath(native, "SPICE", "lsk", "kernel.tls"), "stand-in kernel")
        lib = joinpath(native, "libGRAM.fake")
        write(lib, "stand-in native library")
        paths = Dict("gram-root" => native, "lib" => lib, "spice" => joinpath(native, "SPICE"))
        r = GRecipe.read_recipe(RECIPE)
        delete!(r["provenance"], "native_library_sha256")
        delete!(r["provenance"], "input_lock_file")
        delete!(r["provenance"], "input_lock_sha256")
        r["grid"]["altitude_bounds_km"] = [100.0, 101.0]
        r["grid"]["latitude_bounds_deg"] = [40.0, 41.0]
        r["grid"]["longitude_step_deg"] = 180.0
        recipe = joinpath(tmp, "recipe.toml"); GRecipe.write_toml(recipe, r)
        dry = joinpath(tmp, "dry")
        GRecipe.execute(recipe, dry; package=pkg)
        @test read(joinpath(dry, "recipe.toml")) == read(recipe)
        @test TOML.parsefile(joinpath(dry, "receipt.toml"))["native_process_started"] === false
        @test !isdir(joinpath(dry, "outputs"))
        assert_manifest(dry)
        @test_throws ArgumentError GRecipe.execute(recipe, dry; package=pkg)
        @test_throws ArgumentError GRecipe.execute(recipe, joinpath(tmp, "invalidmode"); mode="automatic", package=pkg)
        @test_throws ArgumentError GRecipe.execute(recipe, joinpath(tmp, "missingpaths"); mode="native", package=pkg)
        pre = joinpath(tmp, "preflight")
        GRecipe.execute(recipe, pre; mode="preflight", paths=paths, package=pkg)
        pre_record = TOML.parsefile(joinpath(pre, "receipt.toml"))
        @test pre_record["status"] == "preflight_passed_not_executed"
        @test pre_record["native_process_started"] === false
        assets = TOML.parsefile(joinpath(pre, "assets_before.toml"))
        @test assets["planet_data"]["table.bin"]["sha256"] == bytes2hex(sha256("selected Mars input"))
        @test length(assets["planet_data"]) == 1
        @test haskey(assets["spice_directory_inventory"], "lsk/kernel.tls")
        lock = Dict("schema_version" => 1, "planet" => "mars", "source" => "synthetic test fixture", "planet_data" => assets["planet_data"], "spice_directory_inventory" => assets["spice_directory_inventory"])
        @test GRecipe.resolved_assets(r, paths; input_lock=lock) == assets
        changed_lock = deepcopy(lock)
        changed_lock["planet_data"]["table.bin"]["sha256"] = repeat("0", 64)
        @test_throws ArgumentError GRecipe.resolved_assets(r, paths; input_lock=changed_lock)
        write(joinpath(native, "Mars", "data", "extra.bin"), "unrecorded input")
        @test_throws ArgumentError GRecipe.resolved_assets(r, paths; input_lock=lock)
        rm(joinpath(native, "Mars", "data", "extra.bin"))
        assert_manifest(pre)
        success = joinpath(tmp, "success")
        lock_path = joinpath(tmp, "test-inputs.toml")
        GRecipe.write_toml(lock_path, lock)
        r["provenance"]["input_lock_file"] = basename(lock_path)
        r["provenance"]["input_lock_sha256"] = bytes2hex(open(sha256, lock_path))
        GRecipe.write_toml(recipe, r)
        withenv("GRAM_LIB" => "wrong inherited value must not be used", "JULIA_PROJECT" => "/does/not/exist/inherited-project", "JULIA_LOAD_PATH" => "/does/not/exist/inherited-load-path") do
            GRecipe.execute(recipe, success; mode="native", paths=paths, package=pkg)
        end
        receipt = TOML.parsefile(joinpath(success, "receipt.toml"))
        @test receipt["status"] == "generated_not_scientifically_accepted"
        @test receipt["assets_unchanged"]
        @test receipt["child_environment"]["project"]["sha256"] == bytes2hex(sha256("[deps]\n"))
        @test receipt["child_environment"]["manifest_status"] == "absent_no_dependency_resolution"
        @test !isfile(joinpath(success, "child_project", "Manifest.toml"))
        @test receipt["input_identity_policy"] == "require_locked_planet_data_and_spice_inventory"
        @test read(joinpath(success, "test-inputs.toml")) == read(lock_path)
        @test receipt["recipe_epoch_utc"] == "2001-11-07T11:51:04.794789Z"
        @test TOML.parsefile(joinpath(success, "product_validation.toml"))["applied_binding_arguments"]["native_model_version"] == "synthetic stand-in; not GRAM"
        @test occursin("without loading GRAM", read(joinpath(success, "builder.log"), String))
        assert_manifest(success)
        for marker in ("FAIL", "WRONG", "MUTATE")
            write(joinpath(native, marker), "trigger")
            output = joinpath(tmp, lowercase(marker))
            @test_throws Exception GRecipe.execute(recipe, output; mode="native", paths=paths, package=pkg)
            @test !ispath(output)
            failures = filter(p -> startswith(p, basename(output) * ".failed-"), readdir(tmp))
            @test length(failures) == 1
            failed = joinpath(tmp, only(failures))
            @test TOML.parsefile(joinpath(failed, "receipt.toml"))["status"] == "failed"
            assert_manifest(failed)
            rm(joinpath(native, marker))
        end
        write(lib, "version https://git-lfs.github.com/spec/v1\noid sha256:fake\nsize 2\n")
        @test_throws ArgumentError GRecipe.resolved_assets(r, paths)
        write(lock_path, "tampered lock")
        @test_throws ArgumentError GRecipe.read_input_lock(r, recipe)
        write(lib, "")
        @test_throws ArgumentError GRecipe.resolved_assets(r, paths)
        write(lib, "different library")
        r["provenance"]["native_library_sha256"] = repeat("0", 64)
        @test_throws ArgumentError GRecipe.resolved_assets(r, paths)
    end
end
