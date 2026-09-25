# The offline loader and interpolation kernel remain owned by GRAMSuite.jl.
# This file adds native-free construction, validation and an explicit boundary policy.

"""
    resolve_gram_grid_file(planet; surrogate_file="", search_roots=String[])

Find an existing `<planet>_surrogate.jls` payload without initializing GRAM.
An explicit file wins; relative files and search roots are relative to `pwd()`.
Otherwise search the caller's directories in order, then the package's known
grid locations (including the adjacent SpaceAGORA `data/GRAM_surrogate`).
Missing files and undownloaded Git LFS pointers are errors. An invalid file at
a higher-priority location is never silently replaced by a lower-priority file.
This function does not download data or validate its scientific provenance.
"""
function resolve_gram_grid_file(
    planet::AbstractString;
    surrogate_file::AbstractString="",
    search_roots=String[]
)::String
    key = _gram_normalize_planet_name(planet)
    explicit = strip(String(surrogate_file))
    if !isempty(explicit)
        return _gram_grid_file_check(abspath(expanduser(explicit)))
    end
    search_roots isa AbstractString && throw(ArgumentError("search_roots must be a collection of directories, not one string."))
    filename = "$(key)_surrogate.jls"
    candidates = String[]
    for root in search_roots
        root isa AbstractString || throw(ArgumentError("Each search_roots entry must be a directory string."))
        isempty(strip(root)) && throw(ArgumentError("A grid search directory cannot be empty."))
        push!(candidates, joinpath(abspath(expanduser(root)), filename))
    end
    package_root = _package_root()
    grid_root = joinpath(package_root, "GRAM Suite 2.0", "simulation", "GRAM", "static_grids")
    append!(candidates, [
        joinpath(package_root, "data", "GRAM_surrogate", filename),
        joinpath(_workspace_root(), "GRAM_surrogate", filename),
        joinpath(package_root, "GRAM Suite 2.0", _gram_planet_dir_name(key), filename),
        joinpath(grid_root, filename),
        joinpath(grid_root, _GRAM_OFFLINE_SURROGATE_PROFILE, "surrogates", filename)
    ])
    unique!(candidates)
    for file in candidates
        (ispath(file) || islink(file)) && return _gram_grid_file_check(file)
    end
    throw(ArgumentError("No downloaded GRAM surrogate grid for '$key'. Supply surrogate_file or search_roots pointing to an authorized grid. Searched: " * join(candidates, ", ")))
end

function _gram_grid_file_check(file::String)::String
    isfile(file) || throw(ArgumentError("GRAM surrogate grid file is missing or is not a regular file: '$file'. Supply a downloaded .jls payload."))
    filesize(file) > 0 || throw(ArgumentError("GRAM surrogate grid file is empty: '$file'."))
    header = open(io -> read(io, min(256, filesize(file))), file, "r")
    pointer_prefix = codeunits("version https://git-lfs.github.com/spec/v1")
    if length(header) >= length(pointer_prefix) && header[1:length(pointer_prefix)] == pointer_prefix
        throw(ArgumentError("GRAM surrogate grid '$file' is a Git LFS pointer, not downloaded data. Retrieve the authorized LFS object or provide the downloaded payload explicitly."))
    end
    return normpath(file)
end

"""
    GRAMGridAtmosphereModel(; planet, surrogate_file="", search_roots=String[],
                           expected_sha256=nothing, above_grid=:error,
                           vacuum_temperature=200.0)

Load a fixed three-dimensional GRAM surrogate snapshot without a native GRAM
installation, SPICE kernels, native handles or native fallback. Reuses the
existing GRAMSuite offline payload loader and trilinear interpolation kernel.
Accepts `surrogate_trilinear` payloads and `full_grid` payloads with format
`spaceagora_gram_static_grid_v1`. If supplied, a surrogate format must be
`spaceagora_gram_surrogate_trilinear_v1`; older surrogates without it remain valid.

Payload axes `alt_km`, `lat_deg`, `lon_deg` become metres/radians; fields are
`density_kgm3`, `temperature_K` and E/N/U winds in m/s. Coordinates must use the
grid's reference surfaces. Existing payloads do not necessarily establish the
altitude/latitude convention or generation epoch: `metadata` retains supplied
top-level metadata (except field arrays/axes), without inventing missing facts.
Grid bounds establish interpolation coverage, not a physical accuracy guarantee.

The payload is loaded once into an owned model snapshot. A new construction
rereads the file, bypassing the legacy path-only cache. `source_sha256` records
the loaded file; optional `expected_sha256` pins an expected artifact. SHA256
integrity alone does not authenticate scientific provenance. Julia-serialized
payloads must come from a trusted source.

Queries outside altitude/latitude coverage throw `DomainError`. Only an explicit
`above_grid=:vacuum` permits zero density above the ceiling, with the specified
positive temperature and zero wind; it never permits lower-bound extrapolation.
Longitude is periodic. Query elapsed time and the wind selector do not change
the stored snapshot: both wind settings return its stored E/N/U wind. This
preserves the existing offline kernel's behavior and does not emulate evolving
native GRAM, stochastic perturbations, or Titan's native low-altitude fallback.

At an exact latitude pole the API returns stored components for the requested
longitude label, without combining them across labels. This does not establish
a unique physical wind vector at the pole; that requires separate validation.

Treat the contained arrays and metadata as read-only; use `deepcopy` for an
independent mutable copy. The native and hybrid constructors are unchanged.
"""
struct GRAMGridAtmosphereModel
    surrogate::GRAMOfflineSurrogate
    source_sha256::String
    above_grid::Symbol
    vacuum_temperature::Float64
    metadata::Dict{String, Any}
end

function GRAMGridAtmosphereModel(;
    planet::AbstractString,
    surrogate_file::AbstractString="",
    search_roots=String[],
    expected_sha256=nothing,
    above_grid::Symbol=:error,
    vacuum_temperature::Real=200.0
)
    key = _gram_normalize_planet_name(planet)
    above_grid in (:error, :vacuum) || throw(ArgumentError("above_grid must be :error or :vacuum; native fallback and extrapolation are unavailable."))
    temperature = Float64(vacuum_temperature)
    isfinite(temperature) && temperature > 0 || throw(ArgumentError("vacuum_temperature must be finite and positive in kelvin."))
    expected = if expected_sha256 === nothing
        nothing
    elseif expected_sha256 isa AbstractString && occursin(r"^[0-9a-fA-F]{64}$", expected_sha256)
        lowercase(String(expected_sha256))
    else
        throw(ArgumentError("expected_sha256 must be a 64-character hexadecimal SHA256 digest."))
    end
    file = resolve_gram_grid_file(key; surrogate_file=surrogate_file, search_roots=search_roots)
    payload, digest = open(file, "r") do io
        before = bytes2hex(SHA.sha256(io))
        expected === nothing || before == expected || throw(ArgumentError("GRAM grid SHA256 mismatch for '$file': expected $expected, got $before."))
        seekstart(io)
        loaded = try
            deserialize(io)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("Cannot read serialized GRAM surrogate payload '$file': $(sprint(showerror, err))"))
        end
        seekstart(io)
        after = bytes2hex(SHA.sha256(io))
        before == after || throw(ArgumentError("GRAM surrogate grid changed while loading '$file'; retry with an immutable payload."))
        loaded, before
    end
    surrogate = try
        _gram_grid_payload_format_check(payload, file)
        _gram_offline_surrogate_from_payload(payload, file, key; allow_full_grid=true)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("Invalid GRAM surrogate payload '$file': $(sprint(showerror, err))"))
    end
    _gram_validate_grid_snapshot(surrogate)
    metadata = Dict{String, Any}(String(k) => deepcopy(v) for (k, v) in payload if String(k) ∉ ("grid", "fields"))
    return GRAMGridAtmosphereModel(surrogate, digest, above_grid, temperature, metadata)
end

# The pure API accepts the builder's two array payloads. Keep the legacy
# native/hybrid loader restricted to surrogate_trilinear through its default.
function _gram_grid_payload_format_check(payload, file::String)
    payload isa Dict || throw(ArgumentError("Expected Dict payload in '$file'."))
    kind = get(payload, "type", "")
    expected = if kind == "full_grid"
        "spaceagora_gram_static_grid_v1"
    elseif kind == "surrogate_trilinear"
        "spaceagora_gram_surrogate_trilinear_v1"
    else
        throw(ArgumentError("Unsupported native-free grid payload type in '$file'."))
    end
    # Older surrogate payloads omitted format; retain that supported case.
    if kind == "full_grid" || haskey(payload, "format")
        get(payload, "format", nothing) == expected || throw(ArgumentError(
            "Payload type '$kind' requires format '$expected' in '$file'."))
    end
    return nothing
end

function _gram_validate_grid_snapshot(grid::GRAMOfflineSurrogate)
    for (name, axis) in (("altitude", grid.alt_nodes_m), ("latitude", grid.lat_nodes_rad), ("longitude", grid.lon_nodes_rad))
        length(axis) >= 2 || throw(ArgumentError("GRAM grid $name axis must contain at least two nodes."))
        all(isfinite, axis) || throw(ArgumentError("GRAM grid $name axis must be finite."))
        all(i -> axis[i] < axis[i + 1], 1:length(axis)-1) || throw(ArgumentError("GRAM grid $name axis must be strictly increasing, without duplicate nodes."))
    end
    first(grid.lat_nodes_rad) >= -pi/2 && last(grid.lat_nodes_rad) <= pi/2 ||
        throw(ArgumentError("GRAM grid latitude nodes must lie in [-90, 90] degrees."))
    first(grid.lon_nodes_rad) >= 0.0 && last(grid.lon_nodes_rad) < 2pi ||
        throw(ArgumentError("GRAM grid longitude nodes must lie in [0, 360) degrees, without a duplicate seam."))
    for (name, field) in (("density", grid.rho), ("temperature", grid.T), ("east wind", grid.wind_e), ("north wind", grid.wind_n), ("up wind", grid.wind_u))
        all(isfinite, field) || throw(ArgumentError("GRAM grid $name values must all be finite."))
    end
    all(x -> x >= 0.0, grid.rho) || throw(ArgumentError("GRAM grid density must be nonnegative in kg/m^3."))
    all(x -> x > 0.0, grid.T) || throw(ArgumentError("GRAM grid temperature must be positive in kelvin."))
    return grid
end

"""
    density_state(model::GRAMGridAtmosphereModel, h, lat, lon, el_time=0.0, wind=true)

Evaluate the stored snapshot at altitude in metres and angles in radians.
Return `(density_kgm3, temperature_K, wind_ENU_ms)`. The returned wind vector is
local east/north/up, as in GRAMSuite's existing surrogate evaluator; callers
needing a different frame must transform it. `el_time` must be finite but is
not applied; `wind` does not select or suppress components of a stored snapshot.
"""
function density_state(
    model::GRAMGridAtmosphereModel,
    h::Real,
    lat::Real,
    lon::Real,
    el_time::Real=0.0,
    wind::Bool=true
)::Tuple{Float64, Float64, SVector{3, Float64}}
    altitude, latitude, longitude, elapsed = Float64(h), Float64(lat), Float64(lon), Float64(el_time)
    all(isfinite, (altitude, latitude, longitude, elapsed)) || throw(DomainError((h, lat, lon, el_time), "Grid atmosphere query coordinates and elapsed time must be finite."))
    grid = model.surrogate
    _gram_axis_segment_checked(grid.lat_nodes_rad, latitude) === nothing &&
        throw(DomainError(lat, "Latitude $(rad2deg(latitude)) degrees is outside the stored GRAM grid: " *
            "$(rad2deg(first(grid.lat_nodes_rad))) to $(rad2deg(last(grid.lat_nodes_rad))) degrees " *
            "($(first(grid.lat_nodes_rad)) to $(last(grid.lat_nodes_rad)) radians). Native fallback is disabled."))
    state = _gram_offline_surrogate_eval(grid, altitude, latitude, longitude)
    state === nothing || return state
    floor_m, ceiling_m = first(grid.alt_nodes_m), last(grid.alt_nodes_m)
    if model.above_grid === :vacuum && altitude > ceiling_m
        return 0.0, model.vacuum_temperature, SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    altitude < floor_m && throw(DomainError(h, "Altitude $(altitude) metres is below the stored GRAM grid floor " *
        "of $(floor_m) metres (grid $(floor_m) to $(ceiling_m) metres). Queries below the grid are not supported, " *
        "and native fallback is disabled."))
    throw(DomainError(h, "Altitude $(altitude) metres is above the stored GRAM grid ceiling of $(ceiling_m) metres " *
        "(grid $(floor_m) to $(ceiling_m) metres). Native fallback is disabled. Where the model is constructed " *
        "directly, above_grid=:vacuum selects zero density above the ceiling."))
end
