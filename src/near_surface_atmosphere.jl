include("mars_near_surface_scalars.jl")

"""
    GRAMNearSurfaceAtmosphereModel(; planet="Mars", surrogate_file, expected_sha256=nothing)

Native-free frozen Mars atmosphere near the surface. It evaluates density, temperature and pressure from 5 m above the
local terrain up to the payload's areoid-height top, from the component fields of a trusted
`spaceagora_mars_near_surface_scalars_v1` payload (for example the preset `mars_global_near_surface_p20_frozen_v1`).
No native GRAM installation, SPICE kernel or GRAM data file is used.

The payload must be given explicitly. It is hashed before and after deserialization: `expected_sha256` pins it, and a
file that changes while loading is rejected. Supply only trusted Julia-serialized payloads.

Query with [`near_surface_state`](@ref) (geodetic latitude and east longitude in degrees, ellipsoid height in metres) for
the full result: density, temperature, pressure, gas constant, regime, first level, areoid and surface heights,
planetocentric latitude, and the status of the surface-layer model used. `density_state` gives the fixed-grid calling
convention (height in metres, angles in radians) and returns a zero wind vector, because the payload stores no winds.

Queries outside the supported domain throw `DomainError` naming the reason: planetocentric latitude beyond the
payload's limit, surface height at or above its volcano limit, clearance below its minimum, areoid height above its top,
or an unavailable component at that position. Query elapsed time is not applied; the atmosphere is a frozen snapshot.
Treat the model's arrays and metadata as read-only; concurrent queries may share one model.
"""
struct GRAMNearSurfaceAtmosphereModel
    core::MarsNearSurfaceScalars.NearSurfaceModel
    source_sha256::String
    metadata::Dict{String, Any}
end

function _gram_near_surface_file_check(file::String)::String
    isfile(file) || throw(ArgumentError("Near-surface payload file is missing or is not a regular file: '$file'. Supply a downloaded .jls payload."))
    filesize(file) > 0 || throw(ArgumentError("Near-surface payload file is empty: '$file'."))
    header = open(io -> read(io, min(256, filesize(file))), file, "r")
    pointer_prefix = codeunits("version https://git-lfs.github.com/spec/v1")
    if length(header) >= length(pointer_prefix) && header[1:length(pointer_prefix)] == pointer_prefix
        throw(ArgumentError("Near-surface payload '$file' is a Git LFS pointer, not downloaded data. Retrieve the payload itself."))
    end
    return file
end

function GRAMNearSurfaceAtmosphereModel(;
    planet::AbstractString="Mars",
    surrogate_file::AbstractString,
    expected_sha256=nothing
)
    _gram_normalize_planet_name(planet) == "mars" || throw(ArgumentError("Near-surface payloads are available for Mars only."))
    expected = if expected_sha256 === nothing
        nothing
    elseif expected_sha256 isa AbstractString && occursin(r"^[0-9a-fA-F]{64}$", expected_sha256)
        lowercase(String(expected_sha256))
    else
        throw(ArgumentError("expected_sha256 must be a 64-character hexadecimal SHA256 digest."))
    end
    isempty(strip(String(surrogate_file))) && throw(ArgumentError("surrogate_file is required: a near-surface payload is never searched for or downloaded."))
    file = _gram_near_surface_file_check(abspath(expanduser(String(surrogate_file))))
    payload, digest = open(file, "r") do io
        before = bytes2hex(SHA.sha256(io))
        expected === nothing || before == expected || throw(ArgumentError("Near-surface payload SHA256 mismatch for '$file': expected $expected, got $before."))
        seekstart(io)
        loaded = try
            deserialize(io)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("Cannot read serialized near-surface payload '$file': $(sprint(showerror, err))"))
        end
        seekstart(io)
        after = bytes2hex(SHA.sha256(io))
        before == after || throw(ArgumentError("Near-surface payload changed while loading '$file'; retry with an immutable payload."))
        loaded, before
    end
    payload isa AbstractDict || throw(ArgumentError("Invalid near-surface payload '$file': expected a Dict."))
    core = try
        MarsNearSurfaceScalars.NearSurfaceModel(payload, digest)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("Invalid near-surface payload '$file': $(sprint(showerror, err))"))
    end
    metadata = Dict{String, Any}(k => deepcopy(v) for (k, v) in core.metadata)
    return GRAMNearSurfaceAtmosphereModel(core, digest, metadata)
end

"""
    near_surface_state(model::GRAMNearSurfaceAtmosphereModel, lat_geodetic_deg, lon_east_deg, h_m)

Evaluate the frozen near-surface atmosphere at geodetic latitude and east longitude in degrees and height above the
reference ellipsoid in metres. Returns a named tuple: `density_kgm3`, `temperature_K`, `pressure_Pa`, `gas_constant`,
`regime` (`:D3`, `:D4` or `:D5`), `first_level_km`, `areoid_height_km`, `surface_height_km`,
`planetocentric_latitude_deg` and `surface_layer_model_status` (`"qualified"`, a provisional status, or `"not used"`
at or above the first level). Throws `DomainError` outside the supported domain.
"""
near_surface_state(model::GRAMNearSurfaceAtmosphereModel, lat_geodetic_deg::Real, lon_east_deg::Real, h_m::Real) =
    MarsNearSurfaceScalars.near_surface_state(model.core, lat_geodetic_deg, lon_east_deg, h_m)

"""
    density_state(model::GRAMNearSurfaceAtmosphereModel, h, lat, lon, el_time=0.0, wind=true)

Fixed-grid calling convention for the near-surface model: height above the reference ellipsoid in metres, geodetic
latitude and east longitude in radians. Returns `(density_kgm3, temperature_K, wind_ENU_ms)`; the payload stores no
winds, so the wind vector is always zero, whatever `wind` selects. `el_time` must be finite but is not applied.
Throws `DomainError` outside the supported domain.
"""
function density_state(
    model::GRAMNearSurfaceAtmosphereModel,
    h::Real,
    lat::Real,
    lon::Real,
    el_time::Real=0.0,
    wind::Bool=true
)::Tuple{Float64, Float64, SVector{3, Float64}}
    altitude, latitude, longitude, elapsed = Float64(h), Float64(lat), Float64(lon), Float64(el_time)
    all(isfinite, (altitude, latitude, longitude, elapsed)) || throw(DomainError((h, lat, lon, el_time), "Near-surface query coordinates and elapsed time must be finite."))
    state = MarsNearSurfaceScalars.near_surface_state(model.core, rad2deg(latitude), rad2deg(longitude), altitude)
    return state.density_kgm3, state.temperature_K, SVector{3, Float64}(0.0, 0.0, 0.0)
end
