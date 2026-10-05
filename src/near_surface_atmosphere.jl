include("mars_near_surface_scalars.jl")
include("mars_near_surface_winds.jl")

"""
    GRAMNearSurfaceAtmosphereModel(; planet="Mars", surrogate_file, expected_sha256=nothing)

Native-free frozen Mars atmosphere near the surface. It evaluates density, temperature and pressure from 5 m above the
local terrain up to the payload's areoid-height top, from the component fields of a trusted payload (for example the
preset `mars_global_near_surface_p20_frozen_v1`). A `spaceagora_mars_near_surface_scalars_v1` payload (versions 1.0.0
and 1.1.0) stores no winds. A `spaceagora_mars_near_surface_v2` payload (version 1.2.0) has the same scalars plus a
wind layer: east, north and vertical winds, clipped at 0.7 times the speed of sound as native does
(`MarsNearSurfaceWinds`). No native GRAM installation, SPICE kernel or GRAM data file is used.

The payload must be given explicitly. It is hashed before and after deserialization: `expected_sha256` pins it, and a
file that changes while loading is rejected. Supply only trusted Julia-serialized payloads.

Query with [`near_surface_state`](@ref) (geodetic latitude and east longitude in degrees, ellipsoid height in metres) for
the scalar result: density, temperature, pressure, gas constant, regime, first level, areoid and surface heights,
planetocentric latitude, and the status of the surface-layer model used. With a wind layer,
[`near_surface_wind_state`](@ref) adds the winds and the speed of sound. `density_state` gives the fixed-grid calling
convention (height in metres, angles in radians). Its wind vector is the stored winds when the payload has a wind
layer, and zero otherwise, because a version 1 payload stores no winds.

Queries outside the supported domain throw `DomainError` naming the reason: planetocentric latitude beyond the
payload's limit, surface height at or above its volcano limit, clearance below its minimum, areoid height above its top,
or an unavailable component at that position (including, with a wind layer, a wind node or a speed of sound that is
not determined there). Query elapsed time is not applied; the atmosphere is a frozen snapshot.
Treat the model's arrays and metadata as read-only; concurrent queries may share one model.
"""
struct GRAMNearSurfaceAtmosphereModel
    core::MarsNearSurfaceScalars.NearSurfaceModel
    source_sha256::String
    metadata::Dict{String, Any}
    winds::Union{Nothing, MarsNearSurfaceWinds.WindLayer}
end
GRAMNearSurfaceAtmosphereModel(core, source_sha256, metadata) = GRAMNearSurfaceAtmosphereModel(core, source_sha256, metadata, nothing)

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
    core, winds = try
        c = MarsNearSurfaceScalars.NearSurfaceModel(payload, digest)
        c, get(payload, "format", "") == MarsNearSurfaceScalars.FORMAT_V2 ? MarsNearSurfaceWinds.WindLayer(payload, c) : nothing
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("Invalid near-surface payload '$file': $(sprint(showerror, err))"))
    end
    metadata = Dict{String, Any}(k => deepcopy(v) for (k, v) in core.metadata)
    return GRAMNearSurfaceAtmosphereModel(core, digest, metadata, winds)
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
    near_surface_winds_available(model::GRAMNearSurfaceAtmosphereModel)

Whether the model's payload has a wind layer (format `spaceagora_mars_near_surface_v2`, version 1.2.0).
"""
near_surface_winds_available(model::GRAMNearSurfaceAtmosphereModel) = model.winds !== nothing

"""
    near_surface_wind_state(model::GRAMNearSurfaceAtmosphereModel, lat_geodetic_deg, lon_east_deg, h_m)

The scalar result of [`near_surface_state`](@ref) and, from the payload's wind layer:
- `wind_east_ms`, `wind_north_ms` and `wind_up_ms`: the winds, with the horizontal components clipped at ±0.7 times
  the speed of sound, as native clips them;
- `unclipped_wind_east_ms` and `unclipped_wind_north_ms`: the same winds before clipping;
- `sound_speed_ms`, and `wind_clipped` (whether either horizontal component was clipped);
- `wind_regime`, from `:D1` (81 km top band) to `:D5` (5 to 30 m clearance);
- `composition_side`: `:dry` at or below 80.0 km areoid height, `:switched` above;
- `in_switch_band`: within 1e-9 km of 80.0 km, where native's sound speed changes composition and its side may differ
  from the model's (the two height computations differ in their last bits).

Throws `ArgumentError` for a payload without winds, and `DomainError` where the scalars, a needed wind node or the speed
of sound is unavailable.
"""
function near_surface_wind_state(model::GRAMNearSurfaceAtmosphereModel, lat_geodetic_deg::Real, lon_east_deg::Real, h_m::Real)
    model.winds === nothing && throw(ArgumentError("This near-surface payload stores no winds; version 1.2.0 (format $(MarsNearSurfaceScalars.FORMAT_V2)) has them."))
    st = MarsNearSurfaceScalars.near_surface_state(model.core, lat_geodetic_deg, lon_east_deg, h_m)
    merge(st, MarsNearSurfaceWinds.winds(model.winds, model.core, st, lon_east_deg))
end

"""
    density_state(model::GRAMNearSurfaceAtmosphereModel, h, lat, lon, el_time=0.0, wind=true)

Fixed-grid calling convention for the near-surface model: height above the reference ellipsoid in metres, geodetic
latitude and east longitude in radians. Returns `(density_kgm3, temperature_K, wind_ENU_ms)`.
- With a wind layer, the wind vector is the clipped east, north and vertical wind of [`near_surface_wind_state`](@ref),
  whatever `wind` selects; as for the stored winds of the upper grid presets, the caller masks them.
- Without one, the wind vector is zero, because the payload stores no winds.

`el_time` must be finite but is not applied. Throws `DomainError` outside the supported domain.
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
    model.winds === nothing && return state.density_kgm3, state.temperature_K, SVector{3, Float64}(0.0, 0.0, 0.0)
    w = MarsNearSurfaceWinds.winds(model.winds, model.core, state, rad2deg(longitude))
    return state.density_kgm3, state.temperature_K, SVector{3, Float64}(w.wind_east_ms, w.wind_north_ms, w.wind_up_ms)
end
