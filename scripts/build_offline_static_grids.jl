#!/usr/bin/env julia

using Dates
using Serialization
using Printf
using Statistics

# Advanced offline generator. Native assets are supplied separately through
# SPACEAGORA_GRAM_ROOT; loading this file does not load the native binding.
const CLI_OPTIONS = Set((
    "alt-max-km",
    "alt-min-km",
    "build-grid",
    "build-surrogate",
    "day",
    "earth-merra2-path",
    "elapsed-time-s",
    "equatorial-radius-km",
    "help",
    "hour",
    "is-planetocentric",
    "lat-max-deg",
    "lat-min-deg",
    "lat-nodes-deg",
    "lib",
    "mars-dust-mode",
    "mars-map-year",
    "mars-min-max",
    "mars-mola-heights",
    "minute",
    "month",
    "nalt",
    "nlat",
    "nlon",
    "out-dir",
    "perturbation-action",
    "planets",
    "polar-radius-km",
    "second",
    "seed",
    "spice",
    "start-time-frame",
    "start-time-scale",
    "strict",
    "summary-file",
    "surrogate-adaptive-alt",
    "surrogate-adaptive-max-step-km",
    "surrogate-adaptive-min-step-km",
    "surrogate-adaptive-pilot-alt-step-km",
    "surrogate-adaptive-pilot-lat-step-deg",
    "surrogate-adaptive-pilot-lon-step-deg",
    "surrogate-alt-bands",
    "surrogate-alt-bands-earth",
    "surrogate-alt-bands-jupiter",
    "surrogate-alt-bands-mars",
    "surrogate-alt-bands-neptune",
    "surrogate-alt-bands-titan",
    "surrogate-alt-bands-uranus",
    "surrogate-alt-bands-venus",
    "surrogate-dir",
    "surrogate-lat-step-deg",
    "surrogate-lon-step-deg",
    "surrogate-only",
    "surrogate-precision",
    "surrogate-source",
    "year",
))
const CLI_BOOLEAN_OPTIONS = Set(("help", "strict", "build-grid", "build-surrogate",
    "surrogate-only", "surrogate-adaptive-alt", "is-planetocentric",
    "perturbation-action", "mars-mola-heights"))

function parse_cli(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || throw(ArgumentError("Unsupported argument '$arg'. Use --key=value."))
        body = arg[3:end]
        pair = split(body, "=", limit=2)
        key = String(pair[1])
        key in CLI_OPTIONS || throw(ArgumentError("Unknown option --$key. Use --help for supported options."))
        haskey(opts, key) && throw(ArgumentError("Duplicate option --$key."))
        if length(pair) == 2
            value = String(pair[2])
        elseif i < length(args) && !startswith(args[i + 1], "--")
            i += 1
            value = args[i]
        elseif key in CLI_BOOLEAN_OPTIONS
            value = "true"
        else
            throw(ArgumentError("--$key requires a value."))
        end
        isempty(strip(value)) && throw(ArgumentError("--$key requires a nonempty value."))
        key in CLI_BOOLEAN_OPTIONS && parse_reference_bool(value, key)
        opts[key] = value
        i += 1
    end
    return opts
end

function print_help(io::IO=stdout)
    println(io, "Advanced frozen-atmosphere grid generator")
    println(io, "Usage: julia --startup-file=no scripts/build_offline_static_grids.jl --key=value ...")
    println(io, "Direct generation requires SPACEAGORA_GRAM_ROOT containing Julia/GRAM.jl, the native library and planetary data.")
    println(io, "--lib and --spice override the library and SPICE paths. No native files are distributed by this script.")
    println(io, "Resampling: --build-grid=false --build-surrogate=true --surrogate-source=grid reads <out-dir>/<planet>_grid.jls.")
    println(io, "Missing or incompatible grid sources are errors; resampling never silently calls native GRAM.")
    println(io, "Use --strict=true to stop on the first planet failure. Full grids store Float64; surrogate precision is explicit.")
    println(io, "Every option accepts --key=value or --key value; only booleans accept a bare --key.")
    println(io, "Supported options:")
    foreach(key -> println(io, "  --", key), sort!(collect(CLI_OPTIONS)))
end

"""Copy and validate explicit latitude samples in degrees, without native state."""
function validate_latitude_nodes(nodes)
    nodes isa AbstractVector && length(nodes) >= 2 || throw(ArgumentError(
        "--lat-nodes-deg requires at least two comma-separated latitude nodes."))
    values = Float64.(nodes)
    all(isfinite, values) && all(diff(values) .> 0.0) &&
        -90.0 <= first(values) < last(values) <= 90.0 || throw(ArgumentError(
        "--lat-nodes-deg must be finite, strictly ascending and within [-90, 90] degrees."))
    return values
end

function parse_latitude_nodes(opts)
    haskey(opts, "lat-nodes-deg") || return nothing
    conflicts = filter(key -> haskey(opts, key), (
        "nlat", "lat-min-deg", "lat-max-deg", "surrogate-lat-step-deg",
        "surrogate-adaptive-pilot-lat-step-deg"))
    isempty(conflicts) || throw(ArgumentError(
        "--lat-nodes-deg cannot be combined with " * join("--" .* collect(conflicts), ", ") * ". Remove the conflicting latitude options."))
    fields = split(opts["lat-nodes-deg"], ','; keepempty=true)
    values = [tryparse(Float64, strip(field)) for field in fields]
    any(isnothing, values) && throw(ArgumentError(
        "--lat-nodes-deg requires a comma-separated list of numeric degrees, with no empty entries."))
    return validate_latitude_nodes(Float64.(values))
end

# Existing checkouts keep native assets under <package>/GRAM Suite 2.0.
# Public native-free installs use an explicit external root for generation.
const LEGACY_GRAM_ROOT = normpath(joinpath(@__DIR__, "..", "GRAM Suite 2.0"))
const GRAM_ROOT = haskey(ENV, "SPACEAGORA_GRAM_ROOT") ?
    abspath(ENV["SPACEAGORA_GRAM_ROOT"]) : LEGACY_GRAM_ROOT
const SUPPORTED_PLANETS = ("earth", "mars", "venus", "jupiter", "uranus", "neptune", "titan")
const cfg_earth_merra2_path = Ref{Union{Nothing, String}}(nothing)
const PLANET_BODIES = Dict(planet => Symbol("BODY_", uppercase(planet)) for planet in SUPPORTED_PLANETS)

function load_native_binding!(root::String)
    isdefined(@__MODULE__, :GRAM) && return nothing
    binding = joinpath(root, "Julia", "GRAM.jl")
    isfile(binding) || throw(ArgumentError(
        "Native GRAM binding not found at $binding. Set SPACEAGORA_GRAM_ROOT to your external GRAM installation."))
    Base.include(@__MODULE__, binding)
    Core.eval(@__MODULE__, :(using .GRAM))
    return nothing
end

# Tuned defaults for smaller/faster surrogates while preserving low-alt fidelity.
const DEFAULT_SURROGATE_BANDS = Dict(
    "earth" => "0:100:1,100:200:2,200:300:4",
    "mars" => "0:100:1,100:200:2,200:300:4",
    "venus" => "0:100:1,100:200:2,200:300:4",
    "titan" => "0:100:1,100:200:2,200:300:4",
    "jupiter" => "0:200:2,200:600:5,600:2000:20,2000:10000:100",
    "uranus" => "0:200:2,200:600:5,600:2000:20,2000:10000:100",
    "neptune" => "0:200:2,200:600:5,600:2000:20,2000:10000:100"
)

parse_bool(s::AbstractString) = parse_reference_bool(s, "boolean option")
parse_int(s::AbstractString) = parse(Int, strip(String(s)))
parse_float(s::AbstractString) = parse(Float64, strip(String(s)))

# These options describe native generation only. Resampling records the request
# but preserves the source atmosphere and never applies these setters.
function parse_reference_bool(raw::AbstractString, option::String)::Bool
    value = lowercase(strip(raw))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError("--$option must be true or false."))
end

function reference_generation_config(opts)
    radius(name) = haskey(opts, name) ? parse_float(opts[name]) : nothing
    equatorial = radius("equatorial-radius-km")
    polar = radius("polar-radius-km")
    (equatorial === nothing) == (polar === nothing) || throw(ArgumentError(
        "--equatorial-radius-km and --polar-radius-km must be supplied together."))
    if equatorial !== nothing
        all(x -> isfinite(x) && x > 0.0, (equatorial, polar)) ||
            throw(ArgumentError("Reference radii must be finite and positive, in kilometres."))
        equatorial >= polar || throw(ArgumentError("Equatorial radius must be at least the polar radius."))
    end
    scale = parse_int(get(opts, "start-time-scale", "1"))
    frame = parse_int(get(opts, "start-time-frame", "1"))
    scale == 1 || throw(ArgumentError("Only --start-time-scale=1 (UTC) is supported by this reference configuration."))
    frame in (0, 1) || throw(ArgumentError("--start-time-frame must be 0 (planet event time) or 1 (Earth-receive time)."))
    map_year = haskey(opts, "mars-map-year") ? parse_int(opts["mars-map-year"]) : nothing
    map_year === nothing || map_year in (0, 1, 2) || throw(ArgumentError("--mars-map-year must be 0, 1 or 2."))
    # The binding accepts Cint, but the verified reference control is the
    # documented ComputeMinMax default 1. Do not infer other native meanings.
    min_max = haskey(opts, "mars-min-max") ? parse_int(opts["mars-min-max"]) : nothing
    min_max === nothing || min_max == 1 || throw(ArgumentError(
        "Only --mars-min-max=1 (the verified ComputeMinMax default) is supported by this reference configuration."))
    dust_mode = lowercase(strip(get(opts, "mars-dust-mode", "legacy-zero")))
    dust_mode in ("legacy-zero", "native-default") || throw(ArgumentError(
        "--mars-dust-mode must be legacy-zero or native-default (leave native dust settings unchanged)."))
    return (
        is_planetocentric=parse_reference_bool(get(opts, "is-planetocentric", "true"), "is-planetocentric"),
        equatorial_radius_km=equatorial, polar_radius_km=polar,
        start_time_scale=scale, start_time_frame=frame,
        perturbation_action=parse_reference_bool(get(opts, "perturbation-action", "false"), "perturbation-action"),
        mars_mola_heights=haskey(opts, "mars-mola-heights") ? parse_reference_bool(opts["mars-mola-heights"], "mars-mola-heights") : nothing,
        mars_map_year=map_year,
        mars_min_max=min_max,
        mars_dust_mode=dust_mode,
    )
end

function require_reference_setter(name::Symbol)
    isdefined(GRAM, name) || throw(ArgumentError("The loaded GRAM binding lacks requested reference setter $name."))
    return nothing
end

function parse_planets(raw::AbstractString)
    txt = lowercase(strip(String(raw)))
    if isempty(txt) || txt in ("all", "*")
        return collect(SUPPORTED_PLANETS)
    end
    planets = String[]
    for token in split(txt, ",")
        p = lowercase(strip(token))
        p in SUPPORTED_PLANETS || error("Unsupported planet '$p'. Supported: $(SUPPORTED_PLANETS)")
        p in planets || push!(planets, p)
    end
    isempty(planets) && error("No valid planets parsed from '$raw'.")
    return planets
end

function parse_surrogate_source(raw::AbstractString)::Symbol
    x = lowercase(strip(String(raw)))
    if x in ("grid", "from-grid")
        return :grid
    elseif x in ("direct", "gram", "from-gram")
        return :direct
    end
    error("Unsupported --surrogate-source='$raw'. Use one of: grid, direct.")
end

function parse_surrogate_precision(raw::AbstractString)::DataType
    x = lowercase(strip(String(raw)))
    if x in ("float32", "f32", "single")
        return Float32
    elseif x in ("float64", "f64", "double")
        return Float64
    end
    error("Unsupported --surrogate-precision='$raw'. Use one of: float32, float64.")
end

function parse_altitude_bands(spec::AbstractString)
    txt = strip(String(spec))
    isempty(txt) && error("Altitude band spec cannot be empty.")
    bands = NTuple{3, Float64}[]
    for token in split(txt, ",")
        t = strip(token)
        isempty(t) && continue
        parts = split(t, ":")
        length(parts) == 3 || error("Invalid altitude band '$t'. Expected a:b:step (km).")
        a = parse_float(parts[1])
        b = parse_float(parts[2])
        step = parse_float(parts[3])
        all(isfinite, (a, b, step)) || throw(ArgumentError("Altitude band values must be finite."))
        step > 0.0 || error("Altitude band step must be > 0, got $step in '$t'.")
        b > a || error("Altitude band upper bound must be > lower bound in '$t'.")
        push!(bands, (a, b, step))
    end
    isempty(bands) && error("No valid altitude bands parsed from '$spec'.")
    return bands
end

function clip_altitude_bands(bands::Vector{NTuple{3, Float64}}, alt_min_km::Float64, alt_max_km::Float64)
    clipped = NTuple{3, Float64}[]
    for (a, b, step) in bands
        lo = max(a, alt_min_km)
        hi = min(b, alt_max_km)
        if hi > lo
            push!(clipped, (lo, hi, step))
        end
    end
    isempty(clipped) && push!(clipped, (alt_min_km, alt_max_km, max(1e-3, alt_max_km - alt_min_km)))
    sort!(clipped, by=x -> x[1])
    return clipped
end

function altitude_nodes_from_bands(
    bands::Vector{NTuple{3, Float64}},
    alt_min_km::Float64,
    alt_max_km::Float64
)
    clipped = clip_altitude_bands(bands, alt_min_km, alt_max_km)
    nodes = Float64[alt_min_km]
    for (a, b, step) in clipped
        x = max(a, alt_min_km)
        if abs(x - nodes[end]) > 1e-9
            push!(nodes, x)
        end
        while x + step < b - 1e-9
            x += step
            abs(x - nodes[end]) > 1e-9 && push!(nodes, x)
        end
        abs(b - nodes[end]) > 1e-9 && push!(nodes, b)
    end
    abs(alt_max_km - nodes[end]) > 1e-9 && push!(nodes, alt_max_km)

    sort!(nodes)
    dedup = Float64[]
    for x in nodes
        if isempty(dedup) || abs(x - dedup[end]) > 1e-9
            push!(dedup, x)
        end
    end
    return dedup
end

function linear_axis_nodes(min_v::Float64, max_v::Float64, step_v::Float64; periodic::Bool=false)
    step_v > 0.0 || error("Axis step must be > 0, got $step_v.")
    nodes = Float64[]
    x = min_v
    while x <= max_v + 1e-9
        push!(nodes, x)
        x += step_v
    end
    if periodic
        if !isempty(nodes) && abs(nodes[end] - max_v) <= 1e-9
            pop!(nodes)
        end
        isempty(nodes) && push!(nodes, min_v)
    else
        if isempty(nodes) || abs(nodes[end] - max_v) > 1e-9
            push!(nodes, max_v)
        end
    end
    return nodes
end

"""Validated sampling latitude bounds in degrees; old configurations keep global coverage."""
function latitude_bounds(cfg)
    nodes = get(cfg, :lat_nodes_deg, nothing)
    if nodes !== nothing
        values = validate_latitude_nodes(nodes)
        return first(values), last(values)
    end
    lower = Float64(get(cfg, :lat_min_deg, -90.0))
    upper = Float64(get(cfg, :lat_max_deg, 90.0))
    isfinite(lower) && isfinite(upper) && -90.0 <= lower < upper <= 90.0 ||
        throw(ArgumentError("Latitude bounds must be finite and satisfy -90 <= --lat-min-deg < --lat-max-deg <= 90 (degrees)."))
    return lower, upper
end

"""Explicit nodes take precedence over legacy programmatic count/step defaults."""
function latitude_axis(cfg; count=nothing, step=nothing)
    nodes = get(cfg, :lat_nodes_deg, nothing)
    nodes === nothing || return validate_latitude_nodes(nodes)
    bounds = latitude_bounds(cfg)
    return count === nothing ? linear_axis_nodes(bounds..., step; periodic=false) :
        collect(range(bounds..., length=count))
end

function latitude_cli_config(opts)
    nodes = parse_latitude_nodes(opts)
    if nodes === nothing
        return (nlat=parse_int(get(opts, "nlat", "73")),
            lat_min_deg=parse_float(get(opts, "lat-min-deg", "-90.0")),
            lat_max_deg=parse_float(get(opts, "lat-max-deg", "90.0")))
    end
    return (nlat=length(nodes), lat_min_deg=first(nodes), lat_max_deg=last(nodes),
        lat_nodes_deg=nodes)
end

function latitude_sampling_metadata(nodes)
    steps = diff(nodes)
    return Dict{String, Any}("lat_mode" => "explicit_nodes", "lat_nodes_count" => length(nodes),
        "lat_min_deg" => first(nodes), "lat_max_deg" => last(nodes),
        "lat_step_min_deg" => minimum(steps), "lat_step_max_deg" => maximum(steps),
        "axis_authority" => "grid.lat_deg; no uniform latitude step applies")
end

function surrogate_latitude_metadata(cfg, nodes, step)
    get(cfg, :lat_nodes_deg, nothing) === nothing && return Dict("lat_step_deg" => step)
    return latitude_sampling_metadata(nodes)
end

function bands_to_string(bands::Vector{NTuple{3, Float64}})
    parts = String[]
    for (a, b, s) in bands
        push!(parts, @sprintf("%.6g:%.6g:%.6g", a, b, s))
    end
    return join(parts, ",")
end

@inline _clamp01(x::Float64) = clamp(x, 0.0, 1.0)

@inline function _safe_quantile(v::Vector{Float64}, p::Float64)::Float64
    isempty(v) && return 0.0
    return quantile(v, _clamp01(p))
end

function bands_from_nodes(nodes::Vector{Float64}; rel_tol::Float64=1e-6)
    length(nodes) >= 2 || error("Need at least two altitude nodes to form bands.")
    b = NTuple{3, Float64}[]
    i0 = 1
    step_ref = nodes[2] - nodes[1]
    for i in 2:(length(nodes) - 1)
        step_i = nodes[i + 1] - nodes[i]
        same_step = abs(step_i - step_ref) <= rel_tol * max(1.0, abs(step_ref))
        if !same_step
            push!(b, (nodes[i0], nodes[i], step_ref))
            i0 = i
            step_ref = step_i
        end
    end
    push!(b, (nodes[i0], nodes[end], step_ref))
    return b
end

@inline function _interp1_clamped(xnodes::Vector{Float64}, ynodes::Vector{Float64}, x::Float64)::Float64
    n = length(xnodes)
    n == length(ynodes) || error("x/ y profile lengths mismatch.")
    n >= 2 || error("Need at least two profile nodes for interpolation.")
    xq = clamp(x, xnodes[1], xnodes[end])
    i0 = clamp(searchsortedlast(xnodes, xq), 1, n - 1)
    i1 = i0 + 1
    x0 = xnodes[i0]
    x1 = xnodes[i1]
    if x1 == x0
        return ynodes[i0]
    end
    w = clamp((xq - x0) / (x1 - x0), 0.0, 1.0)
    return ynodes[i0] + w * (ynodes[i1] - ynodes[i0])
end

function density_profile_from_grid_payload(grid_payload::Dict{String, Any})
    get(grid_payload, "status", "error") == "ok" || error("Grid payload is not usable (status=$(get(grid_payload, "status", "unknown"))).")
    grid = grid_payload["grid"]
    fields = grid_payload["fields"]
    alt_km = Vector{Float64}(grid["alt_km"])
    rho_grid = fields["density_kgm3"]

    nalt = length(alt_km)
    nalt >= 2 || error("Grid payload requires at least two altitude nodes.")
    logrho = Vector{Float64}(undef, nalt)
    @inbounds for ia in 1:nalt
        slice = rho_grid[ia, :, :]
        s = 0.0
        n = length(slice)
        for v in slice
            s += log(max(Float64(v), 1e-30))
        end
        logrho[ia] = s / max(1, n)
    end
    return alt_km, logrho
end

function density_profile_direct(
    root::String,
    planet::String,
    cfg;
    pilot_alt_step_km::Float64,
    pilot_lat_step_deg::Float64,
    pilot_lon_step_deg::Float64
)
    alt_km = linear_axis_nodes(cfg.alt_min_km, cfg.alt_max_km, pilot_alt_step_km; periodic=false)
    lat_deg = latitude_axis(cfg; step=pilot_lat_step_deg)
    lon_deg = linear_axis_nodes(0.0, 360.0, pilot_lon_step_deg; periodic=true)

    @printf("[%s] Building adaptive pilot profile (%d x %d x %d)\n", planet, length(alt_km), length(lat_deg), length(lon_deg))
    logrho = Vector{Float64}(undef, length(alt_km))

    atmos = atmosphere_for_planet(root, planet)
    try
        configure_common!(atmos, cfg)
        configure_planet_specific!(planet, atmos, cfg)

        @inbounds for ia in eachindex(alt_km)
            if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 8)) == 0
                @printf("[%s]   adaptive pilot altitude slice %d / %d\n", planet, ia, length(alt_km))
            end
            h = alt_km[ia]
            s = 0.0
            n = 0
            for lat in lat_deg
                for lon in lon_deg
                    st = sample_state!(atmos, h, lat, lon, cfg.elapsed_time_s, planet; is_planetocentric=get(cfg, :is_planetocentric, true))
                    s += log(max(st.rho, 1e-30))
                    n += 1
                end
            end
            logrho[ia] = s / max(1, n)
        end
    finally
        close!(atmos)
    end

    return alt_km, logrho
end

function adaptive_altitude_nodes(
    alt_profile_km::Vector{Float64},
    logrho_profile::Vector{Float64};
    min_step_km::Float64,
    max_step_km::Float64
)
    n = length(alt_profile_km)
    n == length(logrho_profile) || error("Adaptive profile length mismatch.")
    n >= 2 || error("Adaptive profile requires at least two points.")
    min_step_km > 0.0 || error("Adaptive min step must be > 0.")
    max_step_km >= min_step_km || error("Adaptive max step must be >= min step.")

    grads = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        if i == 1
            dh = max(1e-9, alt_profile_km[2] - alt_profile_km[1])
            grads[i] = abs((logrho_profile[2] - logrho_profile[1]) / dh)
        elseif i == n
            dh = max(1e-9, alt_profile_km[n] - alt_profile_km[n - 1])
            grads[i] = abs((logrho_profile[n] - logrho_profile[n - 1]) / dh)
        else
            dh = max(1e-9, alt_profile_km[i + 1] - alt_profile_km[i - 1])
            grads[i] = abs((logrho_profile[i + 1] - logrho_profile[i - 1]) / dh)
        end
    end

    rho_lo = _safe_quantile(logrho_profile, 0.15)
    rho_hi = _safe_quantile(logrho_profile, 0.85)
    grad_lo = _safe_quantile(grads, 0.15)
    grad_hi = _safe_quantile(grads, 0.85)
    (rho_hi - rho_lo) < 1e-12 && (rho_hi = rho_lo + 1.0)
    (grad_hi - grad_lo) < 1e-12 && (grad_hi = grad_lo + 1.0)

    score = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        s_rho = _clamp01((logrho_profile[i] - rho_lo) / (rho_hi - rho_lo))
        s_grad = _clamp01((grads[i] - grad_lo) / (grad_hi - grad_lo))
        score[i] = max(s_rho, s_grad)
    end

    alt_min_km = alt_profile_km[1]
    alt_max_km = alt_profile_km[end]
    nodes = Float64[alt_min_km]
    while nodes[end] < alt_max_km - 1e-9
        h = nodes[end]
        s = _interp1_clamped(alt_profile_km, score, h)
        dz = max_step_km - s * (max_step_km - min_step_km)
        dz = clamp(dz, min_step_km, max_step_km)
        nxt = min(h + dz, alt_max_km)
        if nxt <= h + 1e-10
            nxt = min(h + min_step_km, alt_max_km)
        end
        push!(nodes, nxt)
    end

    dedup = Float64[]
    for x in nodes
        if isempty(dedup) || abs(x - dedup[end]) > 1e-9
            push!(dedup, x)
        end
    end
    return dedup
end

function resolve_surrogate_bands(
    opts::Dict{String, String},
    planet::String,
    alt_min_km::Float64,
    alt_max_km::Float64
)
    key_specific = "surrogate-alt-bands-$planet"
    spec = if haskey(opts, key_specific)
        opts[key_specific]
    elseif haskey(opts, "surrogate-alt-bands")
        opts["surrogate-alt-bands"]
    else
        get(DEFAULT_SURROGATE_BANDS, planet, @sprintf("%.6g:%.6g:1", alt_min_km, alt_max_km))
    end
    parsed = parse_altitude_bands(spec)
    return clip_altitude_bands(parsed, alt_min_km, alt_max_km)
end

function data_path_for_planet(root::String, planet::String)
    if planet == "earth"
        return joinpath(root, "Earth", "data")
    elseif planet == "mars"
        return joinpath(root, "Mars", "data")
    end
    return nothing
end

function atmosphere_for_planet(root::String, planet::String)
    body = getfield(GRAM, PLANET_BODIES[planet])
    data_path = data_path_for_planet(root, planet)
    if data_path === nothing
        return create_atmosphere(body)
    end
    return create_atmosphere(body; data_path=data_path)
end

function _find_merra2_file(path::AbstractString, month::Integer)
    file = joinpath(path, "All Mean", @sprintf("MERRA2All_%02d.bin", month))
    return isfile(file) ? file : nothing
end

function resolve_earth_merra2_path(root::String, month::Integer; override::Union{Nothing, String}=nothing)
    1 <= month <= 12 || throw(ArgumentError("Earth MERRA2 month must lie in 1:12."))
    explicit = override === nothing ? get(ENV, "EARTH_MERRA2_PATH", nothing) : override
    if explicit !== nothing
        path = normpath(explicit)
        _find_merra2_file(path, month) !== nothing || throw(ArgumentError(
            "The explicit Earth MERRA2 path $path lacks All Mean/MERRA2All_$(@sprintf("%02d", month)).bin; no alternate data directory was selected."))
        return path
    end
    for path in (joinpath(root, "Earth", "data", "MERRA2data"),
                 joinpath(dirname(root), "GRAM_Data", "Earth", "data", "MERRA2data"))
        _find_merra2_file(path, month) !== nothing && return path
    end
    return nothing
end

function configure_common!(atmos, cfg)
    set_start_time!(
        atmos;
        year=cfg.year,
        month=cfg.month,
        day=cfg.day,
        hour=cfg.hour,
        minute=cfg.minute,
        seconds=cfg.second,
        scale=get(cfg, :start_time_scale, 1),
        frame=get(cfg, :start_time_frame, 1)
    )
    set_seed!(atmos, cfg.seed)
    set_perturbation_action!(atmos, get(cfg, :perturbation_action, false))
    set_perturbation_scales!(
        atmos;
        density_scale=0.0,
        ew_wind_scale=0.0,
        ns_wind_scale=0.0,
        vertical_wind_scale=0.0
    )
    return Dict{String, Any}(
        "record_type" => "applied_binding_arguments_not_native_effective_state",
        "unconfigured_native_settings" => "native defaults retained; values and applicability not verified",
        "initial_time" => Dict{String, Any}(String(k) => getproperty(cfg, k) for k in (:year, :month, :day, :hour, :minute, :second)),
        "start_time_scale" => get(cfg, :start_time_scale, 1),
        "start_time_frame" => get(cfg, :start_time_frame, 1),
        "elapsed_time_s" => cfg.elapsed_time_s,
        "seed" => cfg.seed,
        "perturbation_action" => get(cfg, :perturbation_action, false),
        "perturbation_scales" => [0.0, 0.0, 0.0, 0.0],
        "is_planetocentric" => get(cfg, :is_planetocentric, true),
        "position_units" => "kilometres_and_degrees",
        "wind_product" => "stored_mean_ENU",
    )
end

function configure_planet_specific!(planet::String, atmos, cfg=(;))
    applied = Dict{String, Any}()
    equatorial = get(cfg, :equatorial_radius_km, nothing)
    polar = get(cfg, :polar_radius_km, nothing)
    if equatorial !== nothing || polar !== nothing
        planet == "mars" || throw(ArgumentError("Explicit reference radii are currently supported only for Mars."))
        equatorial !== nothing && polar !== nothing || throw(ArgumentError("Both reference radii must be supplied."))
        require_reference_setter(:set_planetary_radii!)
        set_planetary_radii!(atmos; equatorial=equatorial, polar=polar)
        applied["equatorial_radius_km"] = equatorial
        applied["polar_radius_km"] = polar
    end
    if planet == "earth"
        if cfg_earth_merra2_path[] === nothing
            error("Earth MERRA2 data not found. Set --earth-merra2-path=<.../MERRA2data> or EARTH_MERRA2_PATH.")
        end
        earth_set_merra2_path!(atmos, cfg_earth_merra2_path[]::String)
    end
    if planet == "mars"
        dust_mode = get(cfg, :mars_dust_mode, "legacy-zero")
        dust_mode in ("legacy-zero", "native-default") || throw(ArgumentError("Unsupported Mars dust mode: $dust_mode"))
        applied["mars_dust_mode"] = dust_mode
        map_year = get(cfg, :mars_map_year, nothing)
        if map_year !== nothing
            require_reference_setter(:set_map_year!)
            set_map_year!(atmos, map_year)
            applied["mars_map_year"] = map_year
        end
        # set_map_year! leaves the new year pending in the verified native
        # build. set_min_max! pushes its parameters to the Mars sub-models.
        min_max = get(cfg, :mars_min_max, nothing)
        if min_max !== nothing
            min_max isa Integer && !(min_max isa Bool) && min_max == 1 || throw(ArgumentError(
                "Only mars_min_max=1 (the verified ComputeMinMax default) is supported by this reference configuration."))
            require_reference_setter(:set_min_max!)
            applicable(set_min_max!, atmos, min_max) || throw(ArgumentError(
                "The loaded GRAM binding does not support set_min_max!(atmosphere, mars_min_max)."))
            set_min_max!(atmos, min_max)
            applied["mars_min_max"] = min_max
        elseif map_year !== nothing && dust_mode == "native-default"
            @warn "Mars MapYear is requested without a verified parameter push. Use --mars-min-max=1 with --mars-map-year for the Odyssey reference; the requested year may otherwise remain pending in the native sub-models."
        end
        mola = get(cfg, :mars_mola_heights, nothing)
        if mola !== nothing
            require_reference_setter(:set_mola_heights!)
            set_mola_heights!(atmos, mola)
            applied["mars_mola_heights"] = mola
        end
        if dust_mode == "legacy-zero" && isdefined(GRAM, :set_mgcm_dust_levels!)
            set_mgcm_dust_levels!(atmos; constant_level=0.0, min_level=0.0, max_level=0.0)
            applied["mars_mgcm_dust_levels"] = [0.0, 0.0, 0.0]
        end
        if dust_mode == "legacy-zero" && isdefined(GRAM, :set_dust_storm!)
            set_dust_storm!(
                atmos;
                longitude_sun=0.0,
                duration=0.0,
                intensity=0.0,
                max_radius=0.0,
                latitude=0.0,
                longitude=0.0
            )
            applied["mars_dust_storm"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        end
        if dust_mode == "legacy-zero" && isdefined(GRAM, :set_dust_density!)
            set_dust_density!(atmos; nu=0.0, diameter=0.0, density=0.0)
            applied["mars_dust_microphysics"] = [0.0, 0.0, 0.0]
        end
        applied["dust_interpretation"] = dust_mode == "native-default" ?
            "native dust settings left unchanged; effective values and applicability unverified" :
            "legacy off label does not establish a dust-free atmosphere; effective seasonal/TES behavior unverified"
    end
    return applied
end

@inline function sample_state!(
    atmos,
    h_km::Float64,
    lat_deg::Float64,
    lon_deg::Float64,
    elapsed_time_s::Float64,
    planet::String;
    is_planetocentric::Bool=true
)
    set_position!(
        atmos;
        height=h_km,
        latitude=lat_deg,
        longitude=lon_deg,
        elapsed_time=elapsed_time_s,
        is_planetocentric=is_planetocentric
    )
    err = update!(atmos)
    err == 0 || error("GRAM update failed (planet=$planet, alt_km=$h_km, lat_deg=$lat_deg, lon_deg=$lon_deg): $(get_error_message())")

    dyn = get_dynamics_state(atmos)
    wnd = get_winds_state(atmos)
    return (
        rho=Float64(dyn.density),
        temp=Float64(dyn.temperature),
        wind_ew=Float64(wnd.ewWind),
        wind_ns=Float64(wnd.nsWind),
        wind_up=Float64(wnd.verticalWind)
    )
end

function native_model_metadata(planet::String, atmos)
    name = planet == "mars" ? :get_version_string : Symbol(planet, "_get_version_string")
    if !isdefined(GRAM, name) || !applicable(getfield(GRAM, name), atmos)
        return Dict{String, Any}("native_model_version_status" => "getter_not_exposed_by_binding")
    end
    version = strip(replace(String(getfield(GRAM, name)(atmos)), '\0' => ""))
    return Dict{String, Any}(
        "native_model_version_status" => isempty(version) ? "native_getter_returned_empty" : "native_getter",
        "native_model_version" => version,
    )
end

function build_planet_grid(root::String, planet::String, cfg)
    lat_deg = latitude_axis(cfg; count=cfg.nlat)
    @printf("[%s] Building full grid (%d x %d x %d)\n", planet, cfg.nalt, length(lat_deg), cfg.nlon)
    atmos = atmosphere_for_planet(root, planet)
    try
        applied_config = configure_common!(atmos, cfg)
        merge!(applied_config, configure_planet_specific!(planet, atmos, cfg))
        merge!(applied_config, native_model_metadata(planet, atmos))

        alt_km = collect(range(cfg.alt_min_km, cfg.alt_max_km, length=cfg.nalt))
        lon_deg = collect(range(0.0, 360.0, length=cfg.nlon + 1))[1:end-1]
        dims = (cfg.nalt, length(lat_deg), cfg.nlon)

        rho = Array{Float64, 3}(undef, dims...)
        temp = Array{Float64, 3}(undef, dims...)
        wind_ew = Array{Float64, 3}(undef, dims...)
        wind_ns = Array{Float64, 3}(undef, dims...)
        wind_up = Array{Float64, 3}(undef, dims...)

        t0_ns = time_ns()
        @inbounds for ia in eachindex(alt_km)
            if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 10)) == 0
                @printf("[%s]   full-grid altitude slice %d / %d\n", planet, ia, length(alt_km))
            end
            h = alt_km[ia]
            for ilat in eachindex(lat_deg)
                lat = lat_deg[ilat]
                for ilon in eachindex(lon_deg)
                    lon = lon_deg[ilon]
                    s = sample_state!(atmos, h, lat, lon, cfg.elapsed_time_s, planet; is_planetocentric=get(cfg, :is_planetocentric, true))
                    rho[ia, ilat, ilon] = s.rho
                    temp[ia, ilat, ilon] = s.temp
                    wind_ew[ia, ilat, ilon] = s.wind_ew
                    wind_ns[ia, ilat, ilon] = s.wind_ns
                    wind_up[ia, ilat, ilon] = s.wind_up
                end
            end
        end

        elapsed_s = (time_ns() - t0_ns) * 1e-9
        @printf("[%s] Full grid done in %.3f s\n", planet, elapsed_s)

        payload = Dict{String, Any}(
            "status" => "ok",
            "format" => "spaceagora_gram_static_grid_v1",
            "type" => "full_grid",
            "planet" => planet,
            "created_at_utc" => string(now(UTC)),
            "elapsed_s" => elapsed_s,
            "wind_mode" => "on_mean",
            "monte_carlo" => "off",
            "dust" => planet == "mars" && get(cfg, :mars_dust_mode, "legacy-zero") == "native-default" ? "native_default_unmodified" : "off",
            "elapsed_time_s" => cfg.elapsed_time_s,
            "seed" => cfg.seed,
            "initial_time" => deepcopy(applied_config["initial_time"]),
            "generation_config" => applied_config,
            "grid" => Dict(
                "alt_km" => alt_km,
                "lat_deg" => lat_deg,
                "lon_deg" => lon_deg
            ),
            "fields" => Dict(
                "density_kgm3" => rho,
                "temperature_K" => temp,
                "wind_ew_ms" => wind_ew,
                "wind_ns_ms" => wind_ns,
                "wind_up_ms" => wind_up
            )
        )
        if get(cfg, :lat_nodes_deg, nothing) !== nothing
            payload["latitude_sampling"] = latitude_sampling_metadata(lat_deg)
        end
        return payload
    finally
        close!(atmos)
    end
end

@inline function _axis_segment(nodes::Vector{Float64}, x::Float64)
    n = length(nodes)
    n >= 2 || error("Axis must have at least two nodes.")
    xq = clamp(x, nodes[1], nodes[end])
    i0 = clamp(searchsortedlast(nodes, xq), 1, n - 1)
    i1 = i0 + 1
    x0 = nodes[i0]
    x1 = nodes[i1]
    w = x1 == x0 ? 0.0 : (xq - x0) / (x1 - x0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

@inline function _lon_segment(nodes::Vector{Float64}, lon_deg::Float64)
    n = length(nodes)
    n >= 2 || error("Longitude axis must have at least two nodes.")
    period = 360.0
    xq = mod(lon_deg, period)
    xq < 0 && (xq += period)
    i0 = searchsortedlast(nodes, xq)
    i0 = i0 == 0 ? n : clamp(i0, 1, n)
    i1 = i0 == n ? 1 : i0 + 1
    x0 = nodes[i0]
    x1 = i1 == 1 ? nodes[1] + period : nodes[i1]
    xq_adj = i1 == 1 && xq < x0 ? xq + period : xq
    w = x1 == x0 ? 0.0 : (xq_adj - x0) / (x1 - x0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

@inline _lerp(a::Float64, b::Float64, w::Float64) = a + w * (b - a)

@inline function _trilerp(
    c000::Float64, c100::Float64, c010::Float64, c110::Float64,
    c001::Float64, c101::Float64, c011::Float64, c111::Float64,
    wa::Float64, wb::Float64, wc::Float64
)
    c00 = _lerp(c000, c100, wa)
    c10 = _lerp(c010, c110, wa)
    c01 = _lerp(c001, c101, wa)
    c11 = _lerp(c011, c111, wa)
    c0 = _lerp(c00, c10, wb)
    c1 = _lerp(c01, c11, wb)
    return _lerp(c0, c1, wc)
end

function lookup_payload_state(payload::Dict{String, Any}, alt_km::Float64, lat_deg::Float64, lon_deg::Float64)
    get(payload, "status", "error") == "ok" || error("Input payload is not usable (status=$(get(payload, "status", "unknown"))).")
    grid = payload["grid"]
    fields = payload["fields"]
    alt_nodes = grid["alt_km"]::Vector{Float64}
    lat_nodes = grid["lat_deg"]::Vector{Float64}
    lon_nodes = grid["lon_deg"]::Vector{Float64}

    ia0, ia1, wa = _axis_segment(alt_nodes, alt_km)
    ilat0, ilat1, wb = _axis_segment(lat_nodes, lat_deg)
    ilon0, ilon1, wc = _lon_segment(lon_nodes, lon_deg)

    rho_grid = fields["density_kgm3"]::Array{Float64, 3}
    temp_grid = fields["temperature_K"]::Array{Float64, 3}
    w_ew_grid = fields["wind_ew_ms"]::Array{Float64, 3}
    w_ns_grid = fields["wind_ns_ms"]::Array{Float64, 3}
    w_up_grid = fields["wind_up_ms"]::Array{Float64, 3}

    rho = _trilerp(
        rho_grid[ia0, ilat0, ilon0], rho_grid[ia1, ilat0, ilon0],
        rho_grid[ia0, ilat1, ilon0], rho_grid[ia1, ilat1, ilon0],
        rho_grid[ia0, ilat0, ilon1], rho_grid[ia1, ilat0, ilon1],
        rho_grid[ia0, ilat1, ilon1], rho_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    temp = _trilerp(
        temp_grid[ia0, ilat0, ilon0], temp_grid[ia1, ilat0, ilon0],
        temp_grid[ia0, ilat1, ilon0], temp_grid[ia1, ilat1, ilon0],
        temp_grid[ia0, ilat0, ilon1], temp_grid[ia1, ilat0, ilon1],
        temp_grid[ia0, ilat1, ilon1], temp_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    w_ew = _trilerp(
        w_ew_grid[ia0, ilat0, ilon0], w_ew_grid[ia1, ilat0, ilon0],
        w_ew_grid[ia0, ilat1, ilon0], w_ew_grid[ia1, ilat1, ilon0],
        w_ew_grid[ia0, ilat0, ilon1], w_ew_grid[ia1, ilat0, ilon1],
        w_ew_grid[ia0, ilat1, ilon1], w_ew_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    w_ns = _trilerp(
        w_ns_grid[ia0, ilat0, ilon0], w_ns_grid[ia1, ilat0, ilon0],
        w_ns_grid[ia0, ilat1, ilon0], w_ns_grid[ia1, ilat1, ilon0],
        w_ns_grid[ia0, ilat0, ilon1], w_ns_grid[ia1, ilat0, ilon1],
        w_ns_grid[ia0, ilat1, ilon1], w_ns_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    w_up = _trilerp(
        w_up_grid[ia0, ilat0, ilon0], w_up_grid[ia1, ilat0, ilon0],
        w_up_grid[ia0, ilat1, ilon0], w_up_grid[ia1, ilat1, ilon0],
        w_up_grid[ia0, ilat0, ilon1], w_up_grid[ia1, ilat0, ilon1],
        w_up_grid[ia0, ilat1, ilon1], w_up_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )

    return rho, temp, w_ew, w_ns, w_up
end

function _validate_resampling_axis(axis, name::String; latitude::Bool=false, longitude::Bool=false)
    axis isa Vector{Float64} || throw(ArgumentError("$name must be a Vector{Float64}."))
    length(axis) >= 2 || throw(ArgumentError("$name must contain at least two nodes."))
    all(isfinite, axis) || throw(ArgumentError("$name must contain only finite nodes."))
    all(i -> axis[i] < axis[i + 1], 1:length(axis)-1) ||
        throw(ArgumentError("$name must be strictly increasing, without duplicate nodes."))
    if latitude
        first(axis) >= -90.0 && last(axis) <= 90.0 ||
            throw(ArgumentError("$name must lie in [-90, 90] degrees."))
    elseif longitude
        first(axis) >= 0.0 && last(axis) < 360.0 ||
            throw(ArgumentError("$name must lie in [0, 360) degrees, without a duplicate periodic seam."))
    end
    return axis
end

function _validate_resampling_source(payload::Dict{String, Any}, planet::String)
    get(payload, "status", nothing) == "ok" || throw(ArgumentError("Source grid status must be 'ok'."))
    get(payload, "type", nothing) == "full_grid" || throw(ArgumentError("Resampling requires a full_grid source payload."))
    if haskey(payload, "format") && payload["format"] != "spaceagora_gram_static_grid_v1"
        throw(ArgumentError("Unsupported source full-grid format '$(payload["format"])'."))
    end
    get(payload, "planet", nothing) == planet || throw(ArgumentError("Source grid planet does not match requested '$planet'."))
    grid = get(payload, "grid", nothing)
    fields = get(payload, "fields", nothing)
    grid isa AbstractDict || throw(ArgumentError("Source grid axes must be a dictionary."))
    fields isa AbstractDict || throw(ArgumentError("Source grid fields must be a dictionary."))
    alt = _validate_resampling_axis(get(grid, "alt_km", nothing), "Source altitude axis")
    lat = _validate_resampling_axis(get(grid, "lat_deg", nothing), "Source latitude axis"; latitude=true)
    lon = _validate_resampling_axis(get(grid, "lon_deg", nothing), "Source longitude axis"; longitude=true)
    dims = (length(alt), length(lat), length(lon))
    for name in ("density_kgm3", "temperature_K", "wind_ew_ms", "wind_ns_ms", "wind_up_ms")
        field = get(fields, name, nothing)
        field isa Array{Float64, 3} || throw(ArgumentError("Source field '$name' must be an Array{Float64, 3}."))
        size(field) == dims || throw(ArgumentError("Source field '$name' dimensions must match source axes $dims."))
        all(isfinite, field) || throw(ArgumentError("Source field '$name' must contain only finite values."))
    end
    all(x -> x >= 0.0, fields["density_kgm3"]) || throw(ArgumentError("Source density must be nonnegative."))
    all(x -> x > 0.0, fields["temperature_K"]) || throw(ArgumentError("Source temperature must be positive."))
    return grid
end

"""
Resample a stored full grid without regenerating its atmospheric state.

Atmospheric metadata is copied from the source, including any supplied epoch,
elapsed time, seed, forcing and coordinate conventions. Missing provenance stays
missing. Generation settings in `cfg` are not applied to the atmosphere: a warning
reports that policy, and `build_metadata["requested_config"]` records the request.
`build_metadata["source_metadata"]` preserves all source metadata except arrays
and axes, including the original artifact's creation time and build duration.
The legacy top-level `created_at_utc` and `elapsed_s` describe this new artifact.

Altitude and latitude targets must remain within source coverage. Longitude uses
the existing periodic interpolation convention. No native GRAM calls occur here.
Positive density that rounds to zero at the selected storage precision is rejected
with a request to use Float64. Exact zero density and representable subnormals
remain valid; no density floor is imposed.
"""
function build_surrogate_from_grid_payload(
    grid_payload::Dict{String, Any},
    planet::String,
    cfg,
    s_cfg
)
    source_grid = _validate_resampling_source(grid_payload, planet)
    Tnum = s_cfg.precision
    Tnum in (Float32, Float64) || throw(ArgumentError("Resampling precision must be Float32 or Float64."))
    if get(cfg, :lat_nodes_deg, nothing) === nothing
        isfinite(s_cfg.lat_step_deg) && s_cfg.lat_step_deg > 0.0 || throw(ArgumentError("Target latitude step must be finite and positive."))
    end
    isfinite(s_cfg.lon_step_deg) && s_cfg.lon_step_deg > 0.0 || throw(ArgumentError("Target longitude step must be finite and positive."))
    alt_km = s_cfg.alt_nodes === nothing ? altitude_nodes_from_bands(s_cfg.alt_bands, cfg.alt_min_km, cfg.alt_max_km) : Vector{Float64}(s_cfg.alt_nodes)
    lat_deg = latitude_axis(cfg; step=s_cfg.lat_step_deg)
    lon_deg = linear_axis_nodes(0.0, 360.0, s_cfg.lon_step_deg; periodic=true)
    _validate_resampling_axis(alt_km, "Target altitude axis")
    _validate_resampling_axis(lat_deg, "Target latitude axis"; latitude=true)
    _validate_resampling_axis(lon_deg, "Target longitude axis"; longitude=true)
    for (name, target, source) in (("altitude", alt_km, source_grid["alt_km"]), ("latitude", lat_deg, source_grid["lat_deg"]))
        first(target) >= first(source) && last(target) <= last(source) ||
            throw(ArgumentError("Target $name coverage [$(first(target)), $(last(target))] exceeds source coverage [$(first(source)), $(last(source))]; resampling cannot extend the stored atmosphere."))
    end
    dims = (length(alt_km), length(lat_deg), length(lon_deg))

    @warn "Resampling reuses the source atmosphere. Requested generation epoch, elapsed time and seed are not applied; coordinate, radii, time-convention and atmosphere settings are also unapplied. Missing source provenance remains unknown." planet source_elapsed_time_s=get(grid_payload, "elapsed_time_s", nothing) requested_elapsed_time_s=cfg.elapsed_time_s source_seed=get(grid_payload, "seed", nothing) requested_seed=cfg.seed
    @printf("[%s] Building surrogate from grid (%d x %d x %d)\n", planet, dims[1], dims[2], dims[3])

    rho = Array{Tnum, 3}(undef, dims...)
    temp = Array{Tnum, 3}(undef, dims...)
    wind_ew = Array{Tnum, 3}(undef, dims...)
    wind_ns = Array{Tnum, 3}(undef, dims...)
    wind_up = Array{Tnum, 3}(undef, dims...)

    t0_ns = time_ns()
    @inbounds for ia in eachindex(alt_km)
        if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 10)) == 0
            @printf("[%s]   surrogate altitude slice %d / %d\n", planet, ia, length(alt_km))
        end
        h = alt_km[ia]
        for ilat in eachindex(lat_deg)
            lat = lat_deg[ilat]
            for ilon in eachindex(lon_deg)
                lon = lon_deg[ilon]
                r, t, we, wn, wu = lookup_payload_state(grid_payload, h, lat, lon)
                stored_rho = Tnum(r)
                if r > 0.0 && iszero(stored_rho)
                    throw(ArgumentError("Positive density $r kg/m^3 at altitude=$h km, latitude=$lat deg, longitude=$lon deg rounds to zero in $Tnum. Use --surrogate-precision=float64 to preserve positive density."))
                end
                rho[ia, ilat, ilon] = stored_rho
                temp[ia, ilat, ilon] = Tnum(t)
                wind_ew[ia, ilat, ilon] = Tnum(we)
                wind_ns[ia, ilat, ilon] = Tnum(wn)
                wind_up[ia, ilat, ilon] = Tnum(wu)
            end
        end
    end
    elapsed_s = (time_ns() - t0_ns) * 1e-9
    @printf("[%s] Surrogate-from-grid done in %.3f s\n", planet, elapsed_s)
    alt_steps = length(alt_km) >= 2 ? diff(alt_km) : Float64[]
    all(field -> all(isfinite, field), (rho, temp, wind_ew, wind_ns, wind_up)) ||
        throw(ArgumentError("Resampled fields are not finite at the requested precision."))
    all(x -> x > 0, temp) || throw(ArgumentError("Resampled temperature is not positive at the requested precision."))

    source_metadata = Dict{String, Any}(k => deepcopy(v) for (k, v) in grid_payload if k ∉ ("grid", "fields"))
    # Preserve nested coordinate/unit annotations without copying the old axes
    # or field arrays into the new artifact's source-metadata record.
    for (container, data_keys) in (
        ("grid", ("alt_km", "lat_deg", "lon_deg")),
        ("fields", ("density_kgm3", "temperature_K", "wind_ew_ms", "wind_ns_ms", "wind_up_ms")),
    )
        annotations = Dict(k => deepcopy(v) for (k, v) in grid_payload[container] if k ∉ data_keys)
        isempty(annotations) || (source_metadata[container] = annotations)
    end
    payload = deepcopy(source_metadata)
    # Sampling describes the source artifact, not its atmosphere. The complete
    # record remains in build_metadata.source_metadata; target axes are below.
    delete!(payload, "latitude_sampling")
    created_at_utc = string(now(UTC))
    merge!(payload, Dict{String, Any}(
        "status" => "ok",
        "format" => "spaceagora_gram_surrogate_trilinear_v1",
        "type" => "surrogate_trilinear",
        "source" => "grid",
        "planet" => planet,
        "created_at_utc" => created_at_utc,
        "elapsed_s" => elapsed_s,
        "build_metadata" => Dict{String, Any}(
            "operation" => "resample_existing_full_grid",
            "created_at_utc" => created_at_utc,
            "elapsed_s" => elapsed_s,
            "atmosphere_policy" => "source_state_preserved_without_regeneration",
            "requested_config" => Dict{String, Any}(String(k) => deepcopy(getproperty(cfg, k)) for k in propertynames(cfg)),
            "source_metadata" => source_metadata
        ),
        "surrogate" => Dict(
            "precision" => string(Tnum),
            surrogate_latitude_metadata(cfg, lat_deg, s_cfg.lat_step_deg)...,
            "lon_step_deg" => s_cfg.lon_step_deg,
            "alt_mode" => s_cfg.alt_mode,
            "alt_bands_km" => bands_to_string(s_cfg.alt_bands),
            "alt_nodes_count" => length(alt_km),
            "alt_step_min_km" => isempty(alt_steps) ? 0.0 : minimum(alt_steps),
            "alt_step_max_km" => isempty(alt_steps) ? 0.0 : maximum(alt_steps)
        ),
        "grid" => Dict(
            "alt_km" => alt_km,
            "lat_deg" => lat_deg,
            "lon_deg" => lon_deg
        ),
        "fields" => Dict(
            "density_kgm3" => rho,
            "temperature_K" => temp,
            "wind_ew_ms" => wind_ew,
            "wind_ns_ms" => wind_ns,
            "wind_up_ms" => wind_up
        )
    ))
    return payload
end

function build_surrogate_direct(root::String, planet::String, cfg, s_cfg)
    Tnum = s_cfg.precision
    alt_km = s_cfg.alt_nodes === nothing ? altitude_nodes_from_bands(s_cfg.alt_bands, cfg.alt_min_km, cfg.alt_max_km) : Vector{Float64}(s_cfg.alt_nodes)
    lat_deg = latitude_axis(cfg; step=s_cfg.lat_step_deg)
    lon_deg = linear_axis_nodes(0.0, 360.0, s_cfg.lon_step_deg; periodic=true)
    dims = (length(alt_km), length(lat_deg), length(lon_deg))

    @printf("[%s] Building surrogate direct from GRAM (%d x %d x %d)\n", planet, dims[1], dims[2], dims[3])

    rho = Array{Tnum, 3}(undef, dims...)
    temp = Array{Tnum, 3}(undef, dims...)
    wind_ew = Array{Tnum, 3}(undef, dims...)
    wind_ns = Array{Tnum, 3}(undef, dims...)
    wind_up = Array{Tnum, 3}(undef, dims...)

    atmos = atmosphere_for_planet(root, planet)
    try
        applied_config = configure_common!(atmos, cfg)
        merge!(applied_config, configure_planet_specific!(planet, atmos, cfg))
        merge!(applied_config, native_model_metadata(planet, atmos))
        t0_ns = time_ns()
        @inbounds for ia in eachindex(alt_km)
            if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 10)) == 0
                @printf("[%s]   surrogate-direct altitude slice %d / %d\n", planet, ia, length(alt_km))
            end
            h = alt_km[ia]
            for ilat in eachindex(lat_deg)
                lat = lat_deg[ilat]
                for ilon in eachindex(lon_deg)
                    lon = lon_deg[ilon]
                    s = sample_state!(atmos, h, lat, lon, cfg.elapsed_time_s, planet; is_planetocentric=get(cfg, :is_planetocentric, true))
                    rho[ia, ilat, ilon] = Tnum(s.rho)
                    temp[ia, ilat, ilon] = Tnum(s.temp)
                    wind_ew[ia, ilat, ilon] = Tnum(s.wind_ew)
                    wind_ns[ia, ilat, ilon] = Tnum(s.wind_ns)
                    wind_up[ia, ilat, ilon] = Tnum(s.wind_up)
                end
            end
        end
        elapsed_s = (time_ns() - t0_ns) * 1e-9
        @printf("[%s] Surrogate-direct done in %.3f s\n", planet, elapsed_s)
        alt_steps = length(alt_km) >= 2 ? diff(alt_km) : Float64[]

        return Dict{String, Any}(
            "status" => "ok",
            "format" => "spaceagora_gram_surrogate_trilinear_v1",
            "type" => "surrogate_trilinear",
            "source" => "direct",
            "planet" => planet,
            "created_at_utc" => string(now(UTC)),
            "elapsed_s" => elapsed_s,
            "wind_mode" => "on_mean",
            "monte_carlo" => "off",
            "dust" => planet == "mars" && get(cfg, :mars_dust_mode, "legacy-zero") == "native-default" ? "native_default_unmodified" : "off",
            "elapsed_time_s" => cfg.elapsed_time_s,
            "seed" => cfg.seed,
            "initial_time" => deepcopy(applied_config["initial_time"]),
            "generation_config" => applied_config,
            "surrogate" => Dict(
                "precision" => string(Tnum),
                surrogate_latitude_metadata(cfg, lat_deg, s_cfg.lat_step_deg)...,
                "lon_step_deg" => s_cfg.lon_step_deg,
                "alt_mode" => s_cfg.alt_mode,
                "alt_bands_km" => bands_to_string(s_cfg.alt_bands),
                "alt_nodes_count" => length(alt_km),
                "alt_step_min_km" => isempty(alt_steps) ? 0.0 : minimum(alt_steps),
                "alt_step_max_km" => isempty(alt_steps) ? 0.0 : maximum(alt_steps)
            ),
            "grid" => Dict(
                "alt_km" => alt_km,
                "lat_deg" => lat_deg,
                "lon_deg" => lon_deg
            ),
            "fields" => Dict(
                "density_kgm3" => rho,
                "temperature_K" => temp,
                "wind_ew_ms" => wind_ew,
                "wind_ns_ms" => wind_ns,
                "wind_up_ms" => wind_up
            )
        )
    finally
        close!(atmos)
    end
end

function load_payload(path::String)
    open(path, "r") do io
        payload = deserialize(io)
        payload isa Dict || error("Expected Dict payload in '$path'.")
        return Dict{String, Any}(payload)
    end
end

function save_payload(path::String, payload::Dict{String, Any})
    mkpath(dirname(path))
    open(path, "w") do io
        serialize(io, payload)
    end
    return path
end

toml_escape(s::AbstractString) = replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => " ")

function save_summary(path::String, summary_rows::Vector{Dict{String, Any}}, cfg, opts_meta::Dict{String, Any}, planets)
    lat_min_deg, lat_max_deg = latitude_bounds(cfg)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "created_at_utc = \"$(now(UTC))\"")
        println(io, "format = \"spaceagora_gram_offline_products_v2\"")
        println(io, "atmosphere_settings_role = \"requested_generation_settings; resampled products use source payload metadata\"")
        println(io, "wind = \"on\"")
        println(io, "monte_carlo = \"off\"")
        println(io, "dust = \"$(get(cfg, :mars_dust_mode, "legacy-zero") == "native-default" ? "native_default_unmodified_for_mars" : "off")\"")
        println(io, "elapsed_time_s = $(cfg.elapsed_time_s)")
        println(io, "elapsed_time_s_role = \"requested_generation_setting; resampled products retain source state\"")
        println(io, "full_grid_nalt = $(cfg.nalt)")
        println(io, "full_grid_nlat = $(length(latitude_axis(cfg; count=cfg.nlat)))")
        println(io, "full_grid_nlon = $(cfg.nlon)")
        println(io, "alt_min_km = $(cfg.alt_min_km)")
        println(io, "alt_max_km = $(cfg.alt_max_km)")
        println(io, "lat_min_deg = $lat_min_deg")
        println(io, "lat_max_deg = $lat_max_deg")
        println(io, "latitude_bounds_role = \"requested_sampling_domain; product grid axes are authoritative\"")
        if get(cfg, :lat_nodes_deg, nothing) !== nothing
            println(io, "latitude_mode = \"explicit_nodes\"")
            println(io, "lat_nodes_deg = $(latitude_axis(cfg))")
            println(io, "latitude_step_options_role = \"unused; explicit nodes apply to all generation and resampling latitude axes\"")
        end
        println(io, "planets = [\"$(join(planets, "\", \""))\"]")
        println(io, "build_grid = $(opts_meta["build_grid"])")
        println(io, "build_surrogate = $(opts_meta["build_surrogate"])")
        println(io, "surrogate_only = $(opts_meta["surrogate_only"])")
        println(io, "surrogate_source = \"$(opts_meta["surrogate_source"])\"")
        println(io, "surrogate_precision = \"$(opts_meta["surrogate_precision"])\"")
        println(io, "surrogate_lat_step_deg = $(opts_meta["surrogate_lat_step_deg"])")
        println(io, "surrogate_lon_step_deg = $(opts_meta["surrogate_lon_step_deg"])")
        if haskey(opts_meta, "surrogate_adaptive_alt")
            println(io, "surrogate_adaptive_alt = $(opts_meta["surrogate_adaptive_alt"])")
            println(io, "surrogate_adaptive_min_step_km = $(opts_meta["surrogate_adaptive_min_step_km"])")
            println(io, "surrogate_adaptive_max_step_km = $(opts_meta["surrogate_adaptive_max_step_km"])")
            println(io, "surrogate_adaptive_pilot_alt_step_km = $(opts_meta["surrogate_adaptive_pilot_alt_step_km"])")
            println(io, "surrogate_adaptive_pilot_lat_step_deg = $(opts_meta["surrogate_adaptive_pilot_lat_step_deg"])")
            println(io, "surrogate_adaptive_pilot_lon_step_deg = $(opts_meta["surrogate_adaptive_pilot_lon_step_deg"])")
        end
        println(io)

        for row in summary_rows
            println(io, "[[planet]]")
            println(io, "name = \"$(toml_escape(row["planet"]))\"")
            println(io, "status = \"$(toml_escape(row["status"]))\"")
            println(io, "grid_status = \"$(toml_escape(get(row, "grid_status", "skipped")))\"")
            println(io, "surrogate_status = \"$(toml_escape(get(row, "surrogate_status", "skipped")))\"")
            for key in ("native_model_version", "native_model_version_status")
                haskey(row, key) && println(io, "$key = \"$(toml_escape(row[key]))\"")
            end
            if haskey(row, "grid_file")
                println(io, "grid_file = \"$(toml_escape(row["grid_file"]))\"")
            end
            if haskey(row, "grid_elapsed_s")
                println(io, "grid_elapsed_s = $(row["grid_elapsed_s"])")
            end
            if haskey(row, "surrogate_file")
                println(io, "surrogate_file = \"$(toml_escape(row["surrogate_file"]))\"")
            end
            if haskey(row, "surrogate_elapsed_s")
                println(io, "surrogate_elapsed_s = $(row["surrogate_elapsed_s"])")
            end
            if haskey(row, "surrogate_source")
                println(io, "surrogate_source = \"$(toml_escape(row["surrogate_source"]))\"")
            end
            if get(row, "surrogate_source", nothing) == "grid"
                println(io, "surrogate_atmosphere_policy = \"source_snapshot_preserved\"")
                println(io, "surrogate_source_metadata_authority = \"payload.build_metadata.source_metadata\"")
                source_metadata = get(row, "surrogate_source_metadata", Dict{String, Any}())
                for key in ("elapsed_time_s", "seed")
                    value = get(source_metadata, key, nothing)
                    if value isa Union{Int, Float64} && isfinite(value)
                        println(io, "surrogate_source_$key = $value")
                    else
                        state = haskey(source_metadata, key) ? "see_payload" : "unknown"
                        println(io, "surrogate_source_$(key)_status = \"$state\"")
                    end
                end
                epoch_state = haskey(source_metadata, "initial_time") ? "see_payload" : "unknown"
                println(io, "surrogate_source_epoch_status = \"$epoch_state\"")
            end
            if haskey(row, "surrogate_alt_bands_km")
                println(io, "surrogate_alt_bands_km = \"$(toml_escape(row["surrogate_alt_bands_km"]))\"")
            end
            if haskey(row, "surrogate_alt_mode")
                println(io, "surrogate_alt_mode = \"$(toml_escape(row["surrogate_alt_mode"]))\"")
            end
            if haskey(row, "surrogate_alt_nodes")
                println(io, "surrogate_alt_nodes = $(row["surrogate_alt_nodes"])")
            end
            if haskey(row, "surrogate_alt_step_min_km")
                println(io, "surrogate_alt_step_min_km = $(row["surrogate_alt_step_min_km"])")
            end
            if haskey(row, "surrogate_alt_step_max_km")
                println(io, "surrogate_alt_step_max_km = $(row["surrogate_alt_step_max_km"])")
            end
            if haskey(row, "surrogate_alt_profile_source")
                println(io, "surrogate_alt_profile_source = \"$(toml_escape(row["surrogate_alt_profile_source"]))\"")
            end
            if haskey(row, "error")
                println(io, "error = \"$(toml_escape(row["error"]))\"")
            end
            println(io)
        end
    end
    return path
end

function prepare_build(args::Vector{String})
    opts = parse_cli(args)

    planets = parse_planets(get(opts, "planets", "all"))
    out_dir = get(opts, "out-dir", isfile(joinpath(GRAM_ROOT, "Julia", "GRAM.jl")) ?
        joinpath(GRAM_ROOT, "simulation", "GRAM", "static_grids") : joinpath(pwd(), "static_grids"))
    strict = parse_bool(get(opts, "strict", "false"))

    surrogate_only = parse_bool(get(opts, "surrogate-only", "false"))
    build_grid = parse_bool(get(opts, "build-grid", "true"))
    build_surrogate = parse_bool(get(opts, "build-surrogate", "false"))
    if surrogate_only
        build_grid = false
        build_surrogate = true
        haskey(opts, "surrogate-source") || (opts["surrogate-source"] = "direct")
    end

    cfg = (
        nalt=parse_int(get(opts, "nalt", "81")),
        latitude_cli_config(opts)...,
        nlon=parse_int(get(opts, "nlon", "145")),
        alt_min_km=parse_float(get(opts, "alt-min-km", "0.0")),
        alt_max_km=parse_float(get(opts, "alt-max-km", "2000.0")),
        elapsed_time_s=parse_float(get(opts, "elapsed-time-s", "0.0")),
        year=parse_int(get(opts, "year", "2001")),
        month=parse_int(get(opts, "month", "11")),
        day=parse_int(get(opts, "day", "6")),
        hour=parse_int(get(opts, "hour", "19")),
        minute=parse_int(get(opts, "minute", "0")),
        second=parse_float(get(opts, "second", "32.0")),
        seed=parse_int(get(opts, "seed", "1001")),
        reference_generation_config(opts)...
    )

    latitude_bounds(cfg) # Reject invalid bounds before any native initialization.

    s_source = parse_surrogate_source(get(opts, "surrogate-source", "grid"))
    s_precision = parse_surrogate_precision(get(opts, "surrogate-precision", "float32"))
    s_lat_step_deg = parse_float(get(opts, "surrogate-lat-step-deg", "1.75"))
    s_lon_step_deg = parse_float(get(opts, "surrogate-lon-step-deg", "1.75"))
    s_adaptive_alt = parse_bool(get(opts, "surrogate-adaptive-alt", "false"))
    s_adaptive_min_step_km = parse_float(get(opts, "surrogate-adaptive-min-step-km", "0.5"))
    s_adaptive_max_step_km = parse_float(get(opts, "surrogate-adaptive-max-step-km", "6.0"))
    s_adaptive_pilot_alt_step_km = parse_float(get(opts, "surrogate-adaptive-pilot-alt-step-km", "2.0"))
    s_adaptive_pilot_lat_step_deg = parse_float(get(opts, "surrogate-adaptive-pilot-lat-step-deg", "10.0"))
    s_adaptive_pilot_lon_step_deg = parse_float(get(opts, "surrogate-adaptive-pilot-lon-step-deg", "10.0"))
    s_dir = get(opts, "surrogate-dir", joinpath(out_dir, "surrogates"))

    cfg.nalt >= 2 || error("--nalt must be >= 2")
    cfg.nlat >= 2 || error("--nlat must be >= 2")
    cfg.nlon >= 2 || error("--nlon must be >= 2")
    cfg.alt_max_km > cfg.alt_min_km || error("--alt-max-km must be greater than --alt-min-km")
    s_lat_step_deg > 0.0 || error("--surrogate-lat-step-deg must be > 0")
    0.0 < s_lon_step_deg < 360.0 || error("--surrogate-lon-step-deg must be > 0 and < 360")
    s_adaptive_min_step_km > 0.0 || error("--surrogate-adaptive-min-step-km must be > 0")
    s_adaptive_max_step_km >= s_adaptive_min_step_km || error("--surrogate-adaptive-max-step-km must be >= --surrogate-adaptive-min-step-km")
    s_adaptive_pilot_alt_step_km > 0.0 || error("--surrogate-adaptive-pilot-alt-step-km must be > 0")
    s_adaptive_pilot_lat_step_deg > 0.0 || error("--surrogate-adaptive-pilot-lat-step-deg must be > 0")
    0.0 < s_adaptive_pilot_lon_step_deg < 360.0 || error("--surrogate-adaptive-pilot-lon-step-deg must be > 0 and < 360")

    all(isfinite, (cfg.alt_min_km, cfg.alt_max_km, cfg.elapsed_time_s, cfg.second,
        s_lat_step_deg, s_lon_step_deg, s_adaptive_min_step_km, s_adaptive_max_step_km,
        s_adaptive_pilot_alt_step_km, s_adaptive_pilot_lat_step_deg,
        s_adaptive_pilot_lon_step_deg)) || throw(ArgumentError("All time, altitude and spacing values must be finite."))
    Date(cfg.year, cfg.month, cfg.day) # Validate the calendar without native state.
    0 <= cfg.hour <= 23 && 0 <= cfg.minute <= 59 && 0.0 <= cfg.second < 60.0 ||
        throw(ArgumentError("UTC clock fields must satisfy hour 0:23, minute 0:59, second [0,60)."))
    0 <= cfg.seed <= typemax(Int32) || throw(ArgumentError("--seed must fit a nonnegative 32-bit integer."))
    native_required = build_grid || (build_surrogate && s_source == :direct)
    if native_required
        cfg.equatorial_radius_km === nothing || all(==("mars"), planets) ||
            throw(ArgumentError("Explicit reference radii are currently supported only for Mars."))
        if any(startswith(key, "mars-") for key in keys(opts))
            all(==("mars"), planets) || throw(ArgumentError("Mars settings require --planets=mars; use separate recipes for each planet."))
        end
    end
    # Parse selected band options before loading the native binding.
    for planet in planets
        if build_surrogate && !s_adaptive_alt
            resolve_surrogate_bands(opts, planet, cfg.alt_min_km, cfg.alt_max_km)
        end
    end
    libext = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
    libpath = get(opts, "lib", joinpath(GRAM_ROOT, "Build", "lib", "libGRAM.$libext"))
    spice_path = get(opts, "spice", joinpath(GRAM_ROOT, "SPICE"))
    summary_file = get(opts, "summary-file", "grid_build_summary.toml")
    return (; opts, planets, out_dir, strict, surrogate_only, build_grid, build_surrogate, cfg, s_source, s_precision, s_lat_step_deg, s_lon_step_deg, s_adaptive_alt, s_adaptive_min_step_km, s_adaptive_max_step_km, s_adaptive_pilot_alt_step_km, s_adaptive_pilot_lat_step_deg, s_adaptive_pilot_lon_step_deg, s_dir, libpath, spice_path, summary_file, native_required)
end

function execute_build(plan)
    (; opts, planets, out_dir, strict, surrogate_only, build_grid, build_surrogate, cfg, s_source, s_precision, s_lat_step_deg, s_lon_step_deg, s_adaptive_alt, s_adaptive_min_step_km, s_adaptive_max_step_km, s_adaptive_pilot_alt_step_km, s_adaptive_pilot_lat_step_deg, s_adaptive_pilot_lon_step_deg, s_dir, libpath, spice_path, summary_file, native_required) = plan
    if !build_grid && !build_surrogate
        @warn "Both --build-grid=false and --build-surrogate=false. Nothing to do."
        return nothing
    end
    if native_required
        if "earth" in planets
            earth_merra2_override = haskey(opts, "earth-merra2-path") ? normpath(opts["earth-merra2-path"]) : nothing
            earth_merra2_path = resolve_earth_merra2_path(GRAM_ROOT, cfg.month; override=earth_merra2_override)
            if earth_merra2_path === nothing
                @warn "Earth MERRA2 dataset not found. Earth products will fail unless --earth-merra2-path (or EARTH_MERRA2_PATH) is provided."
            else
                @info "Using Earth MERRA2 path." path=earth_merra2_path
            end
            cfg_earth_merra2_path[] = earth_merra2_path
        end
        set_library!(libpath)
        initialize!(spice_path)
    end
    mkpath(out_dir)
    build_surrogate && mkpath(s_dir)

    summary_rows = Dict{String, Any}[]

    for planet in planets
        row = Dict{String, Any}(
            "planet" => planet,
            "status" => "ok",
            "grid_status" => build_grid ? "pending" : "skipped",
            "surrogate_status" => build_surrogate ? "pending" : "skipped"
        )
        grid_payload = nothing
        grid_file_path = joinpath(out_dir, "$(planet)_grid.jls")

        if build_grid
            try
                grid_payload = build_planet_grid(GRAM_ROOT, planet, cfg)
                save_payload(grid_file_path, grid_payload)
                row["grid_status"] = "ok"
                row["grid_file"] = basename(grid_file_path)
                row["grid_elapsed_s"] = grid_payload["elapsed_s"]
                for key in ("native_model_version", "native_model_version_status")
                    haskey(grid_payload["generation_config"], key) && (row[key] = grid_payload["generation_config"][key])
                end
            catch err
                err_msg = sprint(showerror, err)
                row["status"] = "error"
                row["grid_status"] = "error"
                row["error"] = err_msg
                fail_payload = Dict{String, Any}(
                    "status" => "error",
                    "format" => "spaceagora_gram_static_grid_v1",
                    "type" => "full_grid",
                    "planet" => planet,
                    "created_at_utc" => string(now(UTC)),
                    "error" => err_msg
                )
                save_payload(grid_file_path, fail_payload)
                @warn "Failed to build full grid." planet error=err_msg
                strict && rethrow(err)
            end
        end

        if build_surrogate
            surrogate_file_path = joinpath(s_dir, "$(planet)_surrogate.jls")
            try
                surrogate_payload = nothing
                source_used = s_source
                src = nothing
                src_valid = false

                if s_source == :grid
                    src = grid_payload
                    if !(src isa Dict{String, Any} && get(src, "status", "error") == "ok")
                        if isfile(grid_file_path)
                            loaded = load_payload(grid_file_path)
                            if get(loaded, "status", "error") == "ok"
                                src = loaded
                            end
                        end
                    end
                    src_valid = src isa Dict{String, Any} && get(src, "status", "error") == "ok"
                    src_valid || throw(ArgumentError("Grid source unavailable at $grid_file_path. Supply a valid source or explicitly select --surrogate-source=direct."))
                end

                alt_mode = "bands"
                alt_profile_source = "none"
                if s_adaptive_alt
                    profile_alt = Float64[]
                    profile_logrho = Float64[]
                    if src_valid
                        profile_alt, profile_logrho = density_profile_from_grid_payload(src)
                        alt_profile_source = "grid"
                    else
                        profile_alt, profile_logrho = density_profile_direct(
                            GRAM_ROOT,
                            planet,
                            cfg;
                            pilot_alt_step_km=s_adaptive_pilot_alt_step_km,
                            pilot_lat_step_deg=s_adaptive_pilot_lat_step_deg,
                            pilot_lon_step_deg=s_adaptive_pilot_lon_step_deg
                        )
                        alt_profile_source = "direct"
                    end
                    s_alt_nodes = adaptive_altitude_nodes(
                        profile_alt,
                        profile_logrho;
                        min_step_km=s_adaptive_min_step_km,
                        max_step_km=s_adaptive_max_step_km
                    )
                    # The source-derived profile spans the source altitude axis.
                    # Do not silently label that full range as a requested subset.
                    _validate_resampling_axis(s_alt_nodes, "Adaptive altitude axis")
                    first(s_alt_nodes) == cfg.alt_min_km && last(s_alt_nodes) == cfg.alt_max_km ||
                        throw(ArgumentError("Adaptive altitude nodes span [$(first(s_alt_nodes)), $(last(s_alt_nodes))] km instead of requested [$(cfg.alt_min_km), $(cfg.alt_max_km)] km. Source-grid adaptive resampling requires matching source and requested altitude bounds; use fixed altitude bands for a bounded subset."))
                    s_bands = bands_from_nodes(s_alt_nodes)
                    alt_mode = "adaptive_profile"
                else
                    s_bands = resolve_surrogate_bands(opts, planet, cfg.alt_min_km, cfg.alt_max_km)
                    s_alt_nodes = altitude_nodes_from_bands(s_bands, cfg.alt_min_km, cfg.alt_max_km)
                end

                s_cfg = (
                    alt_bands=s_bands,
                    alt_nodes=s_alt_nodes,
                    alt_mode=alt_mode,
                    lat_step_deg=s_lat_step_deg,
                    lon_step_deg=s_lon_step_deg,
                    precision=s_precision
                )

                if s_source == :grid
                    if src_valid
                        surrogate_payload = build_surrogate_from_grid_payload(src, planet, cfg, s_cfg)
                    end
                else
                    surrogate_payload = build_surrogate_direct(GRAM_ROOT, planet, cfg, s_cfg)
                end

                save_payload(surrogate_file_path, surrogate_payload)
                row["surrogate_status"] = "ok"
                row["surrogate_file"] = basename(surrogate_file_path)
                row["surrogate_elapsed_s"] = surrogate_payload["elapsed_s"]
                row["surrogate_source"] = String(source_used)
                if source_used == :direct
                    for key in ("native_model_version", "native_model_version_status")
                        haskey(surrogate_payload["generation_config"], key) && (row[key] = surrogate_payload["generation_config"][key])
                    end
                end
                if source_used == :grid
                    row["surrogate_source_metadata"] = surrogate_payload["build_metadata"]["source_metadata"]
                end
                row["surrogate_alt_bands_km"] = bands_to_string(s_bands)
                row["surrogate_alt_mode"] = alt_mode
                row["surrogate_alt_nodes"] = length(s_alt_nodes)
                if length(s_alt_nodes) >= 2
                    steps = diff(s_alt_nodes)
                    row["surrogate_alt_step_min_km"] = minimum(steps)
                    row["surrogate_alt_step_max_km"] = maximum(steps)
                end
                row["surrogate_alt_profile_source"] = alt_profile_source
            catch err
                err_msg = sprint(showerror, err)
                row["status"] = "error"
                row["surrogate_status"] = "error"
                row["error"] = haskey(row, "error") ? "$(row["error"]) | surrogate: $err_msg" : err_msg
                fail_payload = Dict{String, Any}(
                    "status" => "error",
                    "format" => "spaceagora_gram_surrogate_trilinear_v1",
                    "type" => "surrogate_trilinear",
                    "planet" => planet,
                    "created_at_utc" => string(now(UTC)),
                    "error" => err_msg
                )
                save_payload(surrogate_file_path, fail_payload)
                @warn "Failed to build surrogate product." planet error=err_msg
                strict && rethrow(err)
            end
        end

        push!(summary_rows, row)
    end

    opts_meta = Dict{String, Any}(
        "build_grid" => build_grid,
        "build_surrogate" => build_surrogate,
        "surrogate_only" => surrogate_only,
        "surrogate_source" => String(s_source),
        "surrogate_precision" => string(s_precision),
        "surrogate_lat_step_deg" => s_lat_step_deg,
        "surrogate_lon_step_deg" => s_lon_step_deg,
        "surrogate_adaptive_alt" => s_adaptive_alt,
        "surrogate_adaptive_min_step_km" => s_adaptive_min_step_km,
        "surrogate_adaptive_max_step_km" => s_adaptive_max_step_km,
        "surrogate_adaptive_pilot_alt_step_km" => s_adaptive_pilot_alt_step_km,
        "surrogate_adaptive_pilot_lat_step_deg" => s_adaptive_pilot_lat_step_deg,
        "surrogate_adaptive_pilot_lon_step_deg" => s_adaptive_pilot_lon_step_deg
    )

    summary_path = joinpath(out_dir, summary_file)
    save_summary(summary_path, summary_rows, cfg, opts_meta, planets)

    println("Wrote summary file: $summary_path")
    for row in summary_rows
        grid_status = get(row, "grid_status", "skipped")
        surr_status = get(row, "surrogate_status", "skipped")
        grid_file = get(row, "grid_file", "-")
        surr_file = get(row, "surrogate_file", "-")
        @printf("  %-8s status=%-5s grid=%-8s (%s) surrogate=%-8s (%s)\n",
            row["planet"], row["status"], grid_status, grid_file, surr_status, surr_file)
    end
    all(row -> row["status"] == "ok", summary_rows) || error("One or more grid products failed; see $summary_path.")
    return summary_rows
end

function main(args::Vector{String}=copy(ARGS))
    opts = parse_cli(args)
    if parse_bool(get(opts, "help", "false"))
        print_help()
        return nothing
    end
    plan = prepare_build(args)
    if plan.native_required
        load_native_binding!(GRAM_ROOT)
        # Binding methods were loaded after main entered its compilation world.
        return Base.invokelatest(execute_build, plan)
    end
    return execute_build(plan)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
