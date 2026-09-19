module GridGenerationRecipe

using Dates
using SHA
using TOML
using UUIDs
using Serialization

const PLANETS = ("earth", "mars", "venus", "jupiter", "uranus", "neptune", "titan")
const HASH_PATTERN = r"^[0-9a-f]{64}$"
const FIELDS = ("density_kgm3", "temperature_K", "wind_ew_ms", "wind_ns_ms", "wind_up_ms")

fail(message) = throw(ArgumentError(message))
function keys_exact(table, required, optional=String[]; label="recipe")
    table isa AbstractDict || fail("$label must be a TOML table")
    extra = setdiff(collect(keys(table)), vcat(required, optional))
    missing = setdiff(required, collect(keys(table)))
    isempty(extra) || fail("Unknown $label keys: $(join(sort(extra), ", "))")
    isempty(missing) || fail("Missing $label keys: $(join(sort(missing), ", "))")
end
isnumber(x) = x isa Real && !(x isa Bool) && isfinite(x)
function integer(x, name; lower=0, upper=typemax(Int32))
    x isa Integer && !(x isa Bool) && lower <= x <= upper || fail("$name must be an integer in $lower:$upper")
    Int(x)
end
function bounds(x, name)
    x isa AbstractVector && length(x) == 2 && all(isnumber, x) && x[1] < x[2] || fail("$name must contain two increasing finite numbers")
    Float64.(x)
end
function interval_count(range, step, name)
    isnumber(step) && step > 0 || fail("$name spacing must be finite and positive")
    n = (range[2] - range[1]) / step
    isfinite(n) && 1 <= n <= typemax(Int32) || fail("$name interval count is unsupported")
    isapprox(n, round(n); atol=1e-10, rtol=1e-12) || fail("$name spacing must divide its range exactly")
    round(Int, n)
end
function parse_epoch(epoch)
    epoch isa String || fail("time.epoch_utc must be a quoted UTC string")
    m = match(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}(?:\.\d{1,9})?)Z$", epoch)
    m === nothing && fail("epoch_utc must use YYYY-MM-DDTHH:MM:SS[.fraction]Z with at most nine fractional digits")
    year, month, day, hour, minute = parse.(Int, m.captures[1:5])
    second = parse(Float64, m.captures[6])
    0 <= hour <= 23 && 0 <= minute <= 59 && 0 <= second < 60 || fail("UTC calendar fields are invalid; leap-second epochs are unsupported")
    try
        Date(year, month, day)
    catch
        fail("UTC calendar date is invalid")
    end
    Dict("year" => year, "month" => month, "day" => day, "hour" => hour, "minute" => minute, "second" => second)
end

"""Validate a direct, deterministic frozen generation recipe without loading GRAM."""
function validate_recipe(r::AbstractDict)
    keys_exact(r, ["schema_version", "id", "version", "planet", "time", "coordinates", "atmosphere", "grid", "provenance"], ["mars"])
    integer(r["schema_version"], "schema_version") == 1 || fail("Unsupported recipe schema_version")
    r["id"] isa String && occursin(r"^[a-z][a-z0-9_]*$", r["id"]) || fail("id must be a lowercase identifier")
    r["version"] isa String && occursin(r"^\d+\.\d+\.\d+$", r["version"]) || fail("version must contain three numeric components")
    r["planet"] in PLANETS || fail("Unsupported planet")
    t = r["time"]
    keys_exact(t, ["epoch_utc", "scale", "scale_code", "frame", "frame_code", "elapsed_seconds"]; label="time")
    parse_epoch(t["epoch_utc"])
    t["scale"] == "UTC" && integer(t["scale_code"], "scale_code") == 1 || fail("Only UTC scale_code=1 is supported")
    t["frame"] == "planet_event_time" && integer(t["frame_code"], "frame_code") == 0 || fail("Frozen recipe requires planet_event_time, frame_code=0")
    isnumber(t["elapsed_seconds"]) && t["elapsed_seconds"] == 0 || fail("Use the frozen UTC epoch with elapsed_seconds=0")
    c = r["coordinates"]
    keys_exact(c, ["latitude", "height", "longitude", "radii_mode"], ["equatorial_radius_km", "polar_radius_km"]; label="coordinates")
    c["latitude"] == "geodetic" && c["height"] == "ellipsoidal_km" && c["longitude"] == "east_positive_periodic_degrees" || fail("Only geodetic latitude, ellipsoidal_km height and east-positive periodic degrees are supported")
    if c["radii_mode"] == "explicit"
        r["planet"] == "mars" || fail("The builder currently applies explicit radii only for Mars")
        a, b = get(c, "equatorial_radius_km", nothing), get(c, "polar_radius_km", nothing)
        all(isnumber, (a, b)) && a >= b > 0 || fail("Both radii must be finite and positive, equatorial >= polar")
    elseif c["radii_mode"] == "native-default"
        r["planet"] != "mars" || fail("Mars recipes require explicit radii")
        any(haskey(c, k) for k in ("equatorial_radius_km", "polar_radius_km")) && fail("native-default radii must not claim explicit values")
    else
        fail("radii_mode must be explicit or native-default")
    end
    a = r["atmosphere"]
    keys_exact(a, ["seed", "perturbation_action", "perturbation_scales", "wind_product"]; label="atmosphere")
    integer(a["seed"], "seed")
    a["perturbation_action"] === false || fail("Frozen generation requires perturbation_action=false")
    ps = a["perturbation_scales"]
    ps isa AbstractVector && length(ps) == 4 && all(x -> isnumber(x) && x == 0, ps) || fail("The builder supports four zero perturbation scales")
    a["wind_product"] == "nominal_mean_ENU" || fail("Unsupported wind product")
    if r["planet"] == "mars"
        haskey(r, "mars") || fail("Mars recipes require the mars table")
        m = r["mars"]
        keys_exact(m, ["map_year", "min_max", "mola_heights", "dust_mode"]; label="mars")
        integer(m["map_year"], "mars.map_year"; upper=2)
        integer(m["min_max"], "mars.min_max") == 1 || fail("MapYear requires the verified min_max=1 parameter push")
        m["mola_heights"] === false || fail("Ellipsoidal height requires mola_heights=false")
        m["dust_mode"] in ("native-default", "legacy-zero") || fail("Unsupported Mars dust mode")
    else
        haskey(r, "mars") && fail("mars settings cannot be supplied for another planet")
    end
    g = r["grid"]
    keys_exact(g, ["altitude_bounds_km", "altitude_step_km", "latitude_bounds_deg", "latitude_step_deg", "longitude_bounds_deg", "longitude_step_deg", "longitude_endpoint", "field_type"]; label="grid")
    alt = bounds(g["altitude_bounds_km"], "altitude")
    lat = bounds(g["latitude_bounds_deg"], "latitude")
    -90 <= lat[1] < lat[2] <= 90 || fail("Latitude must stay within [-90,90]")
    lon = bounds(g["longitude_bounds_deg"], "longitude")
    lon == [0.0, 360.0] && g["longitude_endpoint"] == "periodic_without_duplicate" || fail("Longitude must be full [0,360), without a duplicate endpoint")
    interval_count(alt, g["altitude_step_km"], "altitude")
    interval_count(lat, g["latitude_step_deg"], "latitude")
    interval_count(lon, g["longitude_step_deg"], "longitude") >= 2 || fail("Longitude requires at least two nodes")
    g["field_type"] == "Float64" || fail("This recipe schema uses the accepted uniform-grid Float64 generation path")
    p = r["provenance"]
    keys_exact(p, ["native_model_label", "unprobed_defaults", "pole_wind_policy"], ["native_library_sha256", "reference_grid_sha256", "reference_builder_sha256", "reference_status", "input_lock_file", "input_lock_sha256"]; label="provenance")
    for k in ("native_model_label", "pole_wind_policy")
        p[k] isa String && !isempty(strip(p[k])) || fail("provenance.$k must be nonempty")
    end
    p["unprobed_defaults"] isa AbstractVector && !isempty(p["unprobed_defaults"]) && all(x -> x isa String && !isempty(strip(x)), p["unprobed_defaults"]) || fail("Record version-pinned unprobed defaults explicitly")
    for k in ("native_library_sha256", "reference_grid_sha256", "reference_builder_sha256", "input_lock_sha256")
        !haskey(p, k) || (p[k] isa String && occursin(HASH_PATTERN, p[k])) || fail("provenance.$k must be a lowercase SHA256")
    end
    haskey(p, "input_lock_file") == haskey(p, "input_lock_sha256") || fail("input_lock_file and input_lock_sha256 must be supplied together")
    if haskey(p, "input_lock_file")
        name = p["input_lock_file"]
        name isa String && occursin(r"^[a-zA-Z0-9][a-zA-Z0-9_.-]*\.toml$", name) || fail("input_lock_file must be a simple TOML filename alongside the recipe")
        lowercase(name) in ("recipe.toml", "receipt.toml", "invocation.toml", "requested_settings.toml", "sources_before.toml", "sources_after.toml", "assets_before.toml", "assets_after.toml", "product_validation.toml") && fail("input_lock_file conflicts with a retained run filename")
    end
    r
end
read_recipe(path) = validate_recipe(TOML.parsefile(path))

function builder_options(r, output)
    validate_recipe(r)
    g, c, a, t = r["grid"], r["coordinates"], r["atmosphere"], r["time"]
    epoch = parse_epoch(t["epoch_utc"])
    opts = Dict{String, String}(
        "planets" => r["planet"], "out-dir" => output, "summary-file" => joinpath(output, "summary.toml"),
        "build-grid" => "true", "build-surrogate" => "false", "strict" => "true",
        "nalt" => string(interval_count(g["altitude_bounds_km"], g["altitude_step_km"], "altitude") + 1),
        "nlat" => string(interval_count(g["latitude_bounds_deg"], g["latitude_step_deg"], "latitude") + 1),
        "nlon" => string(interval_count(g["longitude_bounds_deg"], g["longitude_step_deg"], "longitude")),
        "alt-min-km" => string(g["altitude_bounds_km"][1]), "alt-max-km" => string(g["altitude_bounds_km"][2]),
        "lat-min-deg" => string(g["latitude_bounds_deg"][1]), "lat-max-deg" => string(g["latitude_bounds_deg"][2]),
        "start-time-scale" => "1", "start-time-frame" => "0", "elapsed-time-s" => "0.0",
        "is-planetocentric" => "false", "perturbation-action" => "false", "seed" => string(a["seed"]))
    for (k, v) in epoch; opts[k] = string(v); end
    if c["radii_mode"] == "explicit"
        opts["equatorial-radius-km"] = string(c["equatorial_radius_km"])
        opts["polar-radius-km"] = string(c["polar_radius_km"])
    end
    if r["planet"] == "mars"
        for (k, v) in r["mars"]; opts["mars-" * replace(k, "_" => "-")] = string(v); end
    end
    opts
end
argument_vector(opts) = ["--$key=$(opts[key])" for key in sort(collect(keys(opts)))]

function file_identity(path)
    isfile(path) || fail("Required input file is missing: $path")
    filesize(path) > 0 || fail("Required input file is empty: $path")
    open(path) do io
        startswith(String(read(io, min(128, filesize(path)))), "version https://git-lfs.github.com/spec/v1") && fail("Git LFS payload has not been downloaded: $path")
    end
    Dict("sha256" => open(sha256, path) |> bytes2hex, "size_bytes" => filesize(path))
end
function inventory_tree(root; julia_only=false)
    isdir(root) || fail("Required input directory is missing: $root")
    files = Dict{String, Any}()
    for (dir, _, names) in walkdir(root; follow_symlinks=true)
        for name in sort(names)
            julia_only && !endswith(name, ".jl") && continue
            full = joinpath(dir, name)
            files[replace(relpath(full, root), '\\' => '/')] = file_identity(full)
        end
    end
    isempty(files) && fail("Required input directory contains no selected files: $root")
    files
end
function git_identity(path)
    result = Dict{String, Any}("revision" => "unavailable", "tracked_worktree_changes" => "unavailable")
    try
        result["revision"] = readchomp(pipeline(`git --no-optional-locks -C $path rev-parse HEAD`; stderr=devnull))
        result["tracked_worktree_changes"] = !isempty(readchomp(pipeline(`git --no-optional-locks -C $path status --porcelain --untracked-files=no`; stderr=devnull)))
    catch
        # Exported source trees are identified by file hashes, without inventing a revision.
    end
    result
end
function source_inventory(package)
    files = Dict{String, Any}()
    for name in ("Project.toml", "scripts/generate_grid.jl", "scripts/grid_generation_recipe.jl", "scripts/build_offline_static_grids.jl", "src/GRAMSuite.jl", "src/grid_atmosphere.jl")
        path = joinpath(package, name)
        if isfile(path); files[name] = file_identity(path)
        elseif startswith(name, "scripts/"); fail("Missing generator source: $path")
        end
    end
    Dict("package_git" => git_identity(package), "files" => files,
        "wrapper_role" => "source identity only; direct generation calls the external binding, not the wrapper runtime")
end
function read_input_lock(r, recipe_path)
    p = r["provenance"]
    haskey(p, "input_lock_file") || return nothing
    path = joinpath(dirname(abspath(recipe_path)), p["input_lock_file"])
    bytes = read(path)
    bytes2hex(sha256(bytes)) == p["input_lock_sha256"] || fail("Input lock does not match recipe SHA256")
    lock = TOML.parse(String(copy(bytes)))
    keys_exact(lock, ["schema_version", "planet", "source", "planet_data", "spice_directory_inventory"]; label="input lock")
    integer(lock["schema_version"], "input lock schema_version") == 1 && lock["planet"] == r["planet"] || fail("Input lock version or planet mismatch")
    for group in ("planet_data", "spice_directory_inventory")
        lock[group] isa AbstractDict && !isempty(lock[group]) || fail("Input lock $group must be a nonempty file inventory")
        for (name, record) in lock[group]
            keys_exact(record, ["sha256", "size_bytes"]; label="input lock file")
            record["sha256"] isa String && occursin(HASH_PATTERN, record["sha256"]) || fail("Invalid locked file SHA256")
            integer(record["size_bytes"], "locked file size"; lower=1, upper=typemax(Int))
            !isabspath(name) && !any(==(".."), split(replace(name, '\\' => '/'), '/')) || fail("Input lock filenames must be relative and contain no parent traversal")
        end
    end
    (; path, bytes, lock)
end
function resolved_assets(r, paths; input_lock=nothing, check_library_pin=true)
    root = abspath(paths["gram-root"])
    library = abspath(paths["lib"])
    spice = abspath(paths["spice"])
    binding = joinpath(root, "Julia")
    isfile(joinpath(binding, "GRAM.jl")) || fail("GRAM root must contain Julia/GRAM.jl")
    record = Dict{String, Any}(
        "binding" => inventory_tree(binding; julia_only=true),
        "native_library" => file_identity(library),
        "spice_directory_inventory" => inventory_tree(spice),
        "inventory_scope" => "Selected directory contents, not a native opened-file/kernel-load trace")
    p = r["provenance"]
    !check_library_pin || !haskey(p, "native_library_sha256") || record["native_library"]["sha256"] == p["native_library_sha256"] || fail("Native library does not match recipe SHA256; a new version needs a separate recipe and validation")
    if r["planet"] in ("earth", "mars")
        data = joinpath(root, uppercasefirst(r["planet"]), "data")
        record["planet_data"] = inventory_tree(data)
    else
        record["planet_data"] = Dict("scope" => "No external planet-data path is passed by this builder; model tables reside in the identified native library")
    end
    if r["planet"] == "earth"
        haskey(paths, "earth-merra2-path") || fail("Earth recipes require an explicit --earth-merra2-path")
        month = parse_epoch(r["time"]["epoch_utc"])["month"]
        file_identity(joinpath(paths["earth-merra2-path"], "All Mean", "MERRA2All_" * lpad(string(month), 2, '0') * ".bin"))
        record["earth_merra2"] = inventory_tree(abspath(paths["earth-merra2-path"]))
    end
    if input_lock !== nothing
        for group in ("planet_data", "spice_directory_inventory")
            record[group] == input_lock[group] || fail("Selected $group differs from the locked recipe: missing, extra or changed files require a separately identified input set")
        end
    end
    record
end
function write_toml(path, data)
    open(path, "w") do io; TOML.print(io, data; sorted=true); end
end
function write_manifest(root)
    rows = String[]
    for (dir, _, names) in walkdir(root)
        for name in names
            rel = replace(relpath(joinpath(dir, name), root), '\\' => '/')
            rel == "SHA256SUMS" && continue
            push!(rows, bytes2hex(open(sha256, joinpath(dir, name))) * "  " * rel)
        end
    end
    write(joinpath(root, "SHA256SUMS"), join(sort(rows), "\n") * "\n")
end

function validate_product(path, r)
    # Only deserialize the builder output from this explicit generation run.
    p = open(deserialize, path)
    p isa AbstractDict && get(p, "status", nothing) == "ok" && get(p, "planet", nothing) == r["planet"] || fail("Generated payload status or planet is incorrect")
    g = r["grid"]
    expected = (collect(range(g["altitude_bounds_km"]...; length=interval_count(g["altitude_bounds_km"], g["altitude_step_km"], "altitude") + 1)),
        collect(range(g["latitude_bounds_deg"]...; length=interval_count(g["latitude_bounds_deg"], g["latitude_step_deg"], "latitude") + 1)),
        collect(range(0.0, 360.0; length=interval_count(g["longitude_bounds_deg"], g["longitude_step_deg"], "longitude") + 1))[1:end-1])
    for (key, axis) in zip(("alt_km", "lat_deg", "lon_deg"), expected)
        p["grid"][key] == axis || fail("Generated $key axis does not match the recipe")
    end
    dims = length.(expected)
    for key in FIELDS
        values = p["fields"][key]
        eltype(values) == Float64 && size(values) == dims && all(isfinite, values) || fail("Generated $key has invalid shape, type or finite values")
    end
    all(>(0), p["fields"]["density_kgm3"]) && all(>(0), p["fields"]["temperature_K"]) || fail("Generated density and temperature must be positive")
    applied = get(p, "generation_config", nothing)
    applied isa AbstractDict || fail("Generated payload is missing applied binding configuration")
    for (key, wanted) in ("start_time_scale" => 1, "start_time_frame" => 0, "elapsed_time_s" => 0.0, "seed" => r["atmosphere"]["seed"], "is_planetocentric" => false, "perturbation_action" => false, "perturbation_scales" => [0.0, 0.0, 0.0, 0.0], "initial_time" => parse_epoch(r["time"]["epoch_utc"]))
        get(applied, key, nothing) == wanted || fail("Generated binding setting $key does not match the recipe")
    end
    if r["planet"] == "mars"
        for (key, wanted) in r["mars"]
            get(applied, "mars_" * key, nothing) == wanted || fail("Generated Mars $key does not match recipe")
        end
        for key in ("equatorial_radius_km", "polar_radius_km")
            get(applied, key, nothing) == r["coordinates"][key] || fail("Generated $key does not match recipe")
        end
    end
    Dict("applied_binding_arguments" => applied,
        "native_effective_state" => "Not inferred from setter arguments; separate native fidelity checks remain required",
        "payload_identity" => file_identity(path))
end

"""Run dry-run, preflight or explicit native generation. Failed runs retain a sibling directory."""
function execute(recipe_path, output; mode="dry-run", paths=Dict{String, String}(), package=normpath(joinpath(@__DIR__, "..")), executor=run)
    mode in ("dry-run", "preflight", "native") || fail("mode must be dry-run, preflight or native")
    recipe_bytes = read(recipe_path)
    r = validate_recipe(TOML.parse(String(copy(recipe_bytes))))
    original_identity = Dict("sha256" => bytes2hex(sha256(recipe_bytes)), "size_bytes" => length(recipe_bytes))
    output = abspath(output)
    ispath(output) && fail("Output already exists; use a new run directory")
    isdir(dirname(output)) || fail("Output parent directory must already exist")
    !islink(output) || fail("Output may not be a symlink")
    selected = Dict{String, String}(String(k) => abspath(v) for (k,v) in paths)
    if mode != "dry-run"
        all(haskey(selected, k) for k in ("gram-root", "lib", "spice")) || fail("preflight/native requires explicit --gram-root, --lib and --spice")
    end
    staging = joinpath(dirname(output), "." * basename(output) * ".staging-" * string(uuid4()))
    mkdir(staging)
    started = string(now(UTC)) * "Z"
    receipt = Dict{String, Any}("mode" => mode, "status" => "preparing", "started_at_utc" => started, "recipe_epoch_utc" => r["time"]["epoch_utc"], "native_process_started" => false, "julia_version" => string(VERSION), "threads" => 1)
    try
        write(joinpath(staging, "recipe.toml"), recipe_bytes)
        file_identity(recipe_path) == original_identity || fail("Recipe changed after validation")
        file_identity(joinpath(staging, "recipe.toml")) == original_identity || fail("Recipe changed while being copied")
        lock_record = read_input_lock(r, recipe_path)
        input_lock = lock_record === nothing ? nothing : lock_record.lock
        if lock_record !== nothing
            write(joinpath(staging, r["provenance"]["input_lock_file"]), lock_record.bytes)
        end
        receipt["input_identity_policy"] = lock_record === nothing ? "capture_observed_inputs_only" : "require_locked_planet_data_and_spice_inventory"
        child_project = joinpath(staging, "child_project")
        mkdir(child_project)
        child_project_file = joinpath(child_project, "Project.toml")
        write(child_project_file, "[deps]\n")
        child_project_identity = file_identity(child_project_file)
        receipt["child_environment"] = Dict("policy" => "explicit isolated stdlib-only project; no dependency resolution", "project" => child_project_identity, "manifest_status" => "absent_no_dependency_resolution")
        original_source = source_inventory(package)
        write_toml(joinpath(staging, "sources_before.toml"), original_source)
        opts = builder_options(r, joinpath(staging, "outputs"))
        for key in ("lib", "spice", "earth-merra2-path")
            haskey(selected, key) && (opts[key] = selected[key])
        end
        normalized = copy(opts)
        normalized["out-dir"] = "outputs"
        normalized["summary-file"] = "outputs/summary.toml"
        for key in ("lib", "spice", "earth-merra2-path")
            haskey(normalized, key) && (normalized[key] = "<" * key * ">")
        end
        overrides = Dict("JULIA_NUM_THREADS" => "1", "OMP_NUM_THREADS" => "1", "OPENBLAS_NUM_THREADS" => "1", "JULIA_LOAD_PATH" => join(("@", "@stdlib"), Sys.iswindows() ? ';' : ':'))
        haskey(selected, "gram-root") && (overrides["SPACEAGORA_GRAM_ROOT"] = selected["gram-root"])
        haskey(selected, "lib") && (overrides["GRAM_LIB"] = selected["lib"])
        write_toml(joinpath(staging, "invocation.toml"), Dict("builder" => "scripts/build_offline_static_grids.jl", "arguments" => argument_vector(normalized), "julia_command_prefix" => Base.julia_cmd().exec, "additional_process_flags" => ["--startup-file=no", "--threads=1", "--project=child_project"], "working_directory" => "this private run directory", "native_paths" => selected, "environment_overrides" => overrides, "environment_record_scope" => "Explicit generator overrides only; no raw environment captured"))
        write_toml(joinpath(staging, "requested_settings.toml"), r)
        if mode == "dry-run"
            receipt["status"] = "prepared_not_executed"
            receipt["assets_checked"] = false
        else
            before = resolved_assets(r, selected; input_lock=input_lock)
            write_toml(joinpath(staging, "assets_before.toml"), before)
            receipt["assets_checked"] = true
            if mode == "native"
                file_identity(recipe_path) == original_identity || fail("Recipe changed before execution")
                source_inventory(package) == original_source || fail("Generator source changed before execution")
                command = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$child_project $(joinpath(package, "scripts", "build_offline_static_grids.jl")) $(argument_vector(opts))`
                command = Cmd(command; dir=staging)
                command = addenv(command, overrides)
                receipt["native_process_started"] = true
                t = time_ns()
                try
                    open(joinpath(staging, "builder.log"), "w") do log
                        executor(pipeline(command; stdout=log, stderr=log))
                    end
                finally
                    receipt["subprocess_elapsed_seconds"] = (time_ns() - t) / 1e9
                    after = resolved_assets(r, selected; check_library_pin=false)
                    write_toml(joinpath(staging, "assets_after.toml"), after)
                    receipt["assets_unchanged"] = before == after
                end
                receipt["assets_unchanged"] || fail("Native input assets changed during generation")
                product = joinpath(staging, "outputs", r["planet"] * "_grid.jls")
                write_toml(joinpath(staging, "product_validation.toml"), validate_product(product, r))
                receipt["status"] = "generated_not_scientifically_accepted"
            else
                receipt["native_process_started"] = false
                receipt["status"] = "preflight_passed_not_executed"
            end
        end
        file_identity(recipe_path) == original_identity || fail("Recipe changed during this run")
        file_identity(child_project_file) == child_project_identity || fail("Child project changed during this run")
        isfile(joinpath(child_project, "Manifest.toml")) && fail("Unexpected dependency resolution created a child Manifest.toml")
        if lock_record !== nothing
            read(lock_record.path) == lock_record.bytes || fail("Input lock changed during this run")
        end
        final_source = source_inventory(package)
        write_toml(joinpath(staging, "sources_after.toml"), final_source)
        final_source == original_source || fail("Generator source changed during this run")
        receipt["finished_at_utc"] = string(now(UTC)) * "Z"
        write_toml(joinpath(staging, "receipt.toml"), receipt)
        write_manifest(staging)
        mv(staging, output; force=false)
    catch err
        receipt["status"] = "failed"
        receipt["failure_type"] = string(typeof(err))
        receipt["finished_at_utc"] = string(now(UTC)) * "Z"
        write_toml(joinpath(staging, "receipt.toml"), receipt)
        write_manifest(staging)
        failed = output * ".failed-" * string(uuid4())
        mv(staging, failed; force=false)
        println(stderr, "Failed generation evidence retained at ", failed)
        rethrow()
    end
    output
end

function main(args)
    if args == ["--help"]
        println("julia --startup-file=no scripts/generate_grid.jl --recipe=FILE --output=NEW_DIRECTORY [--mode=dry-run|preflight|native] [--gram-root=DIR --lib=FILE --spice=DIR --earth-merra2-path=DIR]")
        return
    end
    opts = Dict{String, String}()
    allowed = ("recipe", "output", "mode", "gram-root", "lib", "spice", "earth-merra2-path")
    for arg in args
        startswith(arg, "--") && occursin('=', arg) || fail("Use explicit --key=value arguments")
        key, value = split(arg[3:end], '='; limit=2)
        key in allowed || fail("Unknown argument: --$key")
        !haskey(opts, key) && !isempty(value) || fail("Duplicate or empty argument: --$key")
        opts[key] = value
    end
    all(haskey(opts, k) for k in ("recipe", "output")) || fail("--recipe and --output are required")
    paths = Dict(k => opts[k] for k in ("gram-root", "lib", "spice", "earth-merra2-path") if haskey(opts, k))
    output = execute(opts["recipe"], opts["output"]; mode=get(opts, "mode", "dry-run"), paths=paths)
    println("Generation record: ", output)
end

end # module
