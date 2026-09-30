"""
Evaluator of native-free near-surface Mars scalars (density, temperature, pressure) from frozen component fields.

Use `GRAMSuite.GRAMNearSurfaceAtmosphereModel` to load a payload; this submodule holds the evaluation rule.

A query (geodetic latitude, east longitude, ellipsoid height) is converted to planetocentric latitude and radius on
the reference ellipsoid. The terrain component gives the query's own MOLA surface height zs and areoid radius exactly
as native Mars-GRAM does (2D linear on its MOLA lattice). The first exposed lower-table level
L1 = max(floor(zs + 0.3) + 1, -5) selects the vertical rule, and the component fields at the four corners of the
lattice cell are interpolated bilinearly:
  D3 (z >= L1): T and R linear, ln p linear between the bracketing table levels; rho = p / (R T);
  D4 (zs + 30 m < z < L1): T linear from T30 to T_L1; p = p_L1 exp((L1 - z) / H), H = (T_L1 + T5) / Q_L1;
  D5 (5 m <= clearance <= 30 m): T linear from T5 to T30; same pressure law; R = R_L1 below L1.
Q_L is evaluated from its fitted form in the native 7.5-degree band and 9-degree surface cell; the result reports that
model's qualification status. Outside the supported domain (|phic| above the limit, zs at or above the volcano limit,
clearance below the minimum, areoid height above the top) and where a needed component is unavailable, a DomainError
names the reason; there is no native fallback.
"""
module MarsNearSurfaceScalars

const FORMAT = "spaceagora_mars_near_surface_scalars_v1"
const _ARRAY_KEYS = ("terrain", "level_T_K", "level_R", "level_lnp", "level_source", "surface_T30_K", "surface_T5_K", "q_models")

struct NearSurfaceModel
    zs::Matrix{Float64}; ra::Matrix{Float64}
    tlat0::Float64; tlon0::Float64; tstep::Float64
    lat0::Float64; lstep::Float64; nlat::Int; nlon::Int
    levels::Vector{Float64}
    T::Array{Float64,3}; R::Array{Float64,3}; lnp::Array{Float64,3}
    T30::Matrix{Float64}; T5::Matrix{Float64}
    q::Dict{NTuple{3,Int},NTuple{9,Float64}}
    qstatus::Dict{NTuple{3,Int},String}
    a_km::Float64; b_km::Float64
    zs_refuse_km::Float64; phic_max_deg::Float64; min_clearance_km::Float64; top_km::Float64
    source_sha256::String
    metadata::Dict{String,Any}
end

_payload_error(message) = throw(ArgumentError(message))
_finite_real(x) = x isa Real && !(x isa Bool) && isfinite(x)

# Bilinear terrain is continuous, so all heights between its extrema can generate L1 values,
# including integer levels that no terrain node uses. Clip latitude before finding the extrema
# so terrain outside the advertised support does not require unused levels.
function _validate_first_levels(levels, zs, t, s)
    limit = min(180.0, Float64(s["phic_max_deg"]) + 1e-9)
    limit < 0 && return
    lat0, step = Float64(t["lat0_deg"]), Float64(t["step_deg"])
    lo, hi = Inf, -Inf
    for lat in (-limit, limit)
        u = (lat - lat0) / step
        i = clamp(floor(Int, u), 0, size(zs, 1) - 2); f = u - i
        for j in axes(zs, 2)
            z = (1 - f) * zs[i+1, j] + f * zs[i+2, j]
            lo, hi = min(lo, z), max(hi, z)
        end
    end
    for i in axes(zs, 1)
        -limit <= lat0 + (i - 1) * step <= limit || continue
        for j in axes(zs, 2)
            lo, hi = min(lo, zs[i, j]), max(hi, zs[i, j])
        end
    end
    # Allow for interpolation roundoff at an integer transition, then apply the strict volcano
    # exclusion. A height at the refusal limit itself cannot produce a supported query.
    padding = 16 * eps(max(abs(lo), abs(hi), 1.0))
    lo -= padding
    hi = min(hi + padding, prevfloat(Float64(s["zs_refuse_km"])))
    lo > hi && return
    first_L1, last_L1 = max(floor(lo + 0.3) + 1, -5.0), max(floor(hi + 0.3) + 1, -5.0)
    count(z -> first_L1 <= z <= last_L1 && isinteger(z), levels) == last_L1 - first_L1 + 1 ||
        _payload_error("levels_km must contain every terrain-derived first level from $first_L1 through $last_L1 km exactly.")
end

"""
    NearSurfaceModel(payload::AbstractDict, source_sha256::AbstractString)

Validate a deserialized `spaceagora_mars_near_surface_scalars_v1` payload and build the evaluator's immutable view.
Component fields may contain finite values or NaN sentinels, but not infinities. The level axis must
contain every supported terrain-derived first level, and cover the advertised top altitude.
The arrays are used as stored; treat them as read-only.
"""
function NearSurfaceModel(d::AbstractDict, digest::AbstractString)
    get(d, "format", "") == FORMAT || _payload_error("Unsupported near-surface payload format; expected $FORMAT.")
    lowercase(string(get(d, "planet", ""))) == "mars" || _payload_error("Near-surface payload planet must be Mars.")
    for key in ("terrain", "lattice", "support", "levels_km", "level_T_K", "level_R", "level_lnp", "surface_T30_K",
                "surface_T5_K", "q_models", "radii_km")
        haskey(d, key) || _payload_error("Near-surface payload is missing '$key'.")
    end
    t, g, s, qm = d["terrain"], d["lattice"], d["support"], d["q_models"]
    zs, ra = t["surface_height_km"], t["areoid_radius_km"]
    zs isa Matrix{Float64} && ra isa Matrix{Float64} && size(zs) == size(ra) && all(>=(2), size(zs)) ||
        _payload_error("Terrain fields must be two equal Float64 matrices.")
    all(isfinite, zs) && all(isfinite, ra) || _payload_error("Terrain fields must be finite.")
    all(_finite_real, (t["lat0_deg"], t["lon0_deg"], t["step_deg"], g["lat0_deg"], g["step_deg"])) && t["step_deg"] > 0 && g["step_deg"] > 0 ||
        _payload_error("Terrain and lattice origins and steps must be finite, with positive steps.")
    nlat, nlon = g["nlat"], g["nlon"]
    nlat isa Integer && nlon isa Integer && nlat >= 2 && nlon >= 2 || _payload_error("Lattice sizes must be integers of at least 2.")
    levels = d["levels_km"]
    levels isa Vector{Float64} && length(levels) >= 2 && all(isfinite, levels) && issorted(levels; lt = <=) ||
        _payload_error("levels_km must be finite and strictly increasing.")
    for key in ("level_T_K", "level_R", "level_lnp")
        a = d[key]
        a isa Array{Float64,3} && size(a) == (length(levels), nlat, nlon) || _payload_error("$key must be a Float64 array of size (levels, nlat, nlon).")
        all(x -> !isinf(x), a) || _payload_error("$key must contain finite values or NaN sentinels.")
    end
    for key in ("surface_T30_K", "surface_T5_K")
        a = d[key]
        a isa Matrix{Float64} && size(a) == (nlat, nlon) || _payload_error("$key must be a Float64 matrix of size (nlat, nlon).")
        all(x -> !isinf(x), a) || _payload_error("$key must contain finite values or NaN sentinels.")
    end
    n = length(qm["band"])
    all(length(qm[k]) == n for k in ("cell", "L", "order", "phic_center", "lam_center")) && size(qm["coef"]) == (n, 6) ||
        _payload_error("Surface-layer model tables must have one row per model and six coefficients.")
    haskey(qm, "status") && length(qm["status"]) != n && _payload_error("Surface-layer model status must have one entry per model.")
    radii = d["radii_km"]
    length(radii) == 2 && all(_finite_real, radii) && 0 < radii[2] <= radii[1] || _payload_error("radii_km must be finite, with 0 < polar <= equatorial.")
    all(_finite_real, (s["zs_refuse_km"], s["phic_max_deg"], s["min_clearance_km"], s["top_areoid_km"])) ||
        _payload_error("Support limits must be finite.")
    s["top_areoid_km"] <= last(levels) || _payload_error("top_areoid_km must not exceed the final stored level.")
    _validate_first_levels(levels, zs, t, s)
    q = Dict{NTuple{3,Int},NTuple{9,Float64}}(); qs = Dict{NTuple{3,Int},String}()
    for n in eachindex(qm["band"])
        key = (qm["band"][n], qm["cell"][n], qm["L"][n])
        haskey(q, key) && _payload_error("Duplicate surface-layer model $key.")
        q[key] = (Float64(qm["order"][n]), qm["phic_center"][n], qm["lam_center"][n], qm["coef"][n, :]...)
        qs[key] = haskey(qm, "status") ? String(qm["status"][n]) : "unrecorded"
    end
    all(v -> all(isfinite, v), values(q)) || _payload_error("Surface-layer model coefficients must be finite.")
    meta = Dict{String,Any}(k => v for (k, v) in d if !(k in _ARRAY_KEYS))
    NearSurfaceModel(zs, ra, t["lat0_deg"], t["lon0_deg"], t["step_deg"], g["lat0_deg"], g["step_deg"], nlat, nlon,
        levels, d["level_T_K"], d["level_R"], d["level_lnp"], d["surface_T30_K"], d["surface_T5_K"], q, qs,
        radii[1], radii[2], s["zs_refuse_km"], s["phic_max_deg"], s["min_clearance_km"], s["top_areoid_km"],
        String(digest), meta)
end

function radial_position(a, b, lat_deg, h_km)
    θ = deg2rad(lat_deg); s, c = sincos(θ); e2 = 1 - (b / a)^2
    n = a / sqrt(1 - e2 * s * s)
    x = (n + h_km) * c; z = (n * (1 - e2) + h_km) * s
    hypot(x, z), rad2deg(atan(z, x))
end

function terrain(m::NearSurfaceModel, F::Matrix{Float64}, phic, lam)
    u = (phic - m.tlat0) / m.tstep
    i = clamp(floor(Int, u), 0, size(F, 1) - 2); fu = u - i
    x = mod(lam - m.tlon0, 360.0) / m.tstep
    # mod rounds up to exactly 360.0 for longitudes within about 3e-14 degrees below the lattice origin, so the column
    # index is wrapped; the value there is the origin column's, which is the continuous limit.
    j = mod(floor(Int, x), size(F, 2)); fv = x - floor(x); j2 = mod(j + 1, size(F, 2))
    (1 - fu) * (1 - fv) * F[i+1, j+1] + fu * (1 - fv) * F[i+2, j+1] + (1 - fu) * fv * F[i+1, j2+1] + fu * fv * F[i+2, j2+1]
end

wrapdeg(d) = mod(d + 180.0, 360.0) - 180.0

function near_surface_state(m::NearSurfaceModel, lat_geodetic_deg::Real, lon_east_deg::Real, h_m::Real)
    all(isfinite, (lat_geodetic_deg, lon_east_deg, h_m)) && -90 <= lat_geodetic_deg <= 90 || throw(DomainError((lat_geodetic_deg, lon_east_deg, h_m), "Finite geodetic position required"))
    r, phic = radial_position(m.a_km, m.b_km, Float64(lat_geodetic_deg), Float64(h_m) / 1000)
    abs(phic) <= m.phic_max_deg + 1e-9 || throw(DomainError(phic, "Planetocentric latitude beyond $(m.phic_max_deg) degrees is not supported"))
    lam = mod(Float64(lon_east_deg), 360.0)
    zs = terrain(m, m.zs, phic, lam); z = r - terrain(m, m.ra, phic, lam)
    zs < m.zs_refuse_km || throw(DomainError(zs, "Surface height $(round(zs; digits=3)) km: volcano flanks (>= $(m.zs_refuse_km) km) are not supported"))
    z - zs >= m.min_clearance_km - 1e-9 || throw(DomainError(z - zs, "Below 5 m above the surface is not supported"))
    z <= m.top_km + 1e-9 || throw(DomainError(z, "Above the $(m.top_km) km areoid-height top of the near-surface product"))
    L1 = max(floor(zs + 0.3) + 1, -5.0)
    ui = (phic - m.lat0) / m.lstep; i = clamp(floor(Int, ui), 0, m.nlat - 2); u = ui - i
    vj = lam / m.lstep; j = mod(floor(Int, vj), m.nlon); v = vj - floor(vj); j2 = mod(j + 1, m.nlon)
    function bil(A::AbstractArray, a)
        w = ((1 - u) * (1 - v), u * (1 - v), (1 - u) * v, u * v)
        vals = (A[a, i+1, j+1], A[a, i+2, j+1], A[a, i+1, j2+1], A[a, i+2, j2+1])
        any(k -> w[k] != 0 && isnan(vals[k]), 1:4) && throw(DomainError((phic, lam, m.levels[a]), "Near-surface component unavailable here (unreconstructed buried level)"))
        sum(w[k] == 0 ? 0.0 : w[k] * vals[k] for k in 1:4)
    end
    function bil2(A::Matrix{Float64})
        w = ((1 - u) * (1 - v), u * (1 - v), (1 - u) * v, u * v)
        vals = (A[i+1, j+1], A[i+2, j+1], A[i+1, j2+1], A[i+2, j2+1])
        any(k -> w[k] != 0 && isnan(vals[k]), 1:4) && throw(DomainError((phic, lam), "Near-surface surface field unavailable here"))
        sum(w[k] == 0 ? 0.0 : w[k] * vals[k] for k in 1:4)
    end
    level(a) = (T = bil(m.T, a), R = bil(m.R, a), lnp = bil(m.lnp, a))
    if z >= L1
        a = searchsortedlast(m.levels, z); a = clamp(a, 1, length(m.levels))
        lo = level(a)
        if a == length(m.levels) || z == m.levels[a]
            T, R, lnp = lo.T, lo.R, lo.lnp
        else
            hi = level(a + 1); f = (z - m.levels[a]) / (m.levels[a+1] - m.levels[a])
            T = lo.T + (hi.T - lo.T) * f; R = lo.R + (hi.R - lo.R) * f; lnp = lo.lnp + (hi.lnp - lo.lnp) * f
        end
        p = exp(lnp)
        return (density_kgm3 = p / (R * T), temperature_K = T, pressure_Pa = p, gas_constant = R, regime = :D3, first_level_km = L1,
                areoid_height_km = z, surface_height_km = zs, planetocentric_latitude_deg = phic, surface_layer_model_status = "not used")
    end
    lv = level(searchsortedfirst(m.levels, L1))
    T30, T5 = bil2(m.T30), bil2(m.T5)
    band, cell = floor(Int, phic / 7.5), mod(floor(Int, lam / 9.0), 40)
    key = (band, cell, Int(L1))
    # A query on a native band or surface-cell edge belongs to both sides; Q_L is continuous there. The alternative is
    # the neighbour across the nearest edge: below it when the query is on or just above the edge, above it when the
    # query is just below.
    if !haskey(m.q, key)
        alts = NTuple{3,Int}[]
        eb = round(Int, phic / 7.5)
        abs(phic / 7.5 - eb) < 1e-9 && push!(alts, (band == eb ? eb - 1 : eb, cell, Int(L1)))
        ec = round(Int, lam / 9.0)
        abs(lam / 9.0 - ec) < 1e-9 && push!(alts, (band, mod(floor(Int, lam / 9.0) == ec ? ec - 1 : ec, 40), Int(L1)))
        k2 = findfirst(k -> haskey(m.q, k), alts)
        k2 === nothing && throw(DomainError((phic, lam, L1), "Surface-layer coefficient unavailable here"))
        key = alts[k2]
    end
    qd = m.q[key]
    x, y = (phic - qd[2]) / 3.75, wrapdeg(lam - qd[3]) / 4.5
    Q = qd[4] + qd[5] * x + qd[6] * x * x + qd[7] * y + qd[8] * x * y + qd[9] * x * x * y
    H = (lv.T + T5) / Q
    z30 = zs + 0.03
    if z > z30
        T = T30 + (lv.T - T30) * (z - z30) / (L1 - z30); regime = :D4
    else
        T = T5 + (T30 - T5) * (z - (zs + 0.005)) / 0.025; regime = :D5
    end
    p = exp(lv.lnp) * exp((L1 - z) / H)
    (density_kgm3 = p / (lv.R * T), temperature_K = T, pressure_Pa = p, gas_constant = lv.R, regime = regime, first_level_km = L1,
     areoid_height_km = z, surface_height_km = zs, planetocentric_latitude_deg = phic, surface_layer_model_status = m.qstatus[key])
end
end # module MarsNearSurfaceScalars
