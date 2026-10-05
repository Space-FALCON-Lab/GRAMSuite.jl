"""
Wind layer of the native-free near-surface Mars atmosphere: east, north and vertical winds, and the sound speed that
limits the horizontal winds, from the stored fields of a `spaceagora_mars_near_surface_v2` payload. Density,
temperature, pressure and gas constant come unchanged from `MarsNearSurfaceScalars`; this module only adds winds.

The query's planetocentric latitude, east longitude, areoid height z and surface height zs are the scalar evaluator's.

**Winds.**
- **Height rule:** five regimes over clearance Δ = z - zs and the first exposed level L1:
  - D5, Δ <= 30 m: the stored 5 m and 30 m surface fields;
  - D4, below L1: the 30 m field and level L1;
  - D3, below 75 km: the bracketing lower-table levels;
  - D2, below 80 + ho km: the 75 km level and the first MTGCM level (80 + ho);
  - D1, above that: the two MTGCM levels (80 + ho and 85 + ho).
- **Horizontal interpolation:** each endpoint field is interpolated linearly in latitude and longitude.
  - The east wind's latitude axis has two-sided knots at 78.75 N and 82.5 N, where native's pole factor jumps.
  - The surface east wind is split at the 0 E seam.
  - Terms with zero weight are skipped, so an unavailable node that the query does not need never contaminates it.
- **Slope and vertical terms:** native's slope-wind term (`MarsWindEndpoints.slope_wind`, unchanged) adds them to the
  background. The terrain slopes come from ±0.25 degree differences of the published terrain, as the fields were
  generated.

**Sound speed.** c = sqrt(γ P / ρ), with γ from a reference K(T, P) (K = γ / (γ - 1)) on native's 50 K and pressure-decade
cells plus the stored offset δ.
- At or below the 80.0 km composition switch, water vapour at native's 20 % relative humidity is mixed back in.
- Above the switch, a composition term at the query's own molecular weight R_u / R is added, with no water.
- In D2 only the endpoint on the query's side of 80.0 km is used.
- The side is chosen by the query's own areoid height. Native's height arithmetic differs in the last bits, so within
  `SWITCH_BAND_KM` of 80.0 km the side may differ from native's.

**Clipping.** The horizontal winds are clipped to ±0.7c. The vertical wind is not clipped.

**Refusals.** A needed node that is not stored, or a sound speed the reference or the composition does not determine,
is refused with a `DomainError` naming it. A wind is never returned as zero for lack of data.
"""
module MarsNearSurfaceWinds

import ..MarsWindEndpoints
import ..MarsNearSurfaceScalars

const FORMAT = MarsNearSurfaceScalars.FORMAT_V2
const ARRAY_KEYS = MarsNearSurfaceScalars.WIND_ARRAY_KEYS
const CLIP_FRACTION = 0.7
const SWITCH_KM = 80.0
const SWITCH_BAND_KM = 1e-9
# The sound-speed reference and its water and composition rules (the accepted W0 composition revision, version 1.1).
const TK = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0]            # temperature levels, K
const PK = [1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0, 1e4, 1e5]            # pressure levels, Pa
const TW = [100.0, 200.0, 300.0, 400.0]                               # NIST-JANAF H2O(g) temperatures, K
const KW = [33.299, 33.349, 33.596, 34.262] ./ 8.314462618            # water vapour's K = Cp / R there
const RH = 0.2                                                        # native's relative humidity
const MD, MW = 43.47, 18.0153                                         # dry-mix and water molecular weights
const RW = MD / MW
const RU = 8314.46261815324                                           # J/(kmol K): M = RU / R

_payload_error(message) = throw(ArgumentError(message))

"A latitude axis: knot latitudes (a two-sided knot appears twice, lower side first) and their sides (-1, +1, or 0)."
struct Axis
    lat::Vector{Float64}
    side::Vector{Int}
end

"The stored fields of a version 2 payload, validated against its scalar part."
struct WindLayer
    U::Array{Float64,3}; V::Array{Float64,3}               # (level, knot, longitude), pre-clip m/s
    MU::Array{Float64,3}; MV::Array{Float64,3}             # (80 + ho or 85 + ho, MTGCM row, longitude)
    SU::Array{Float64,3}; SV::Array{Float64,3}             # (5 m or 30 m, knot, longitude); SU has 241 seam-split longitudes
    O::Array{Float64,3}; OS::Array{Float64,3}              # sound offsets (scalar level, S knot, longitude); (5 m or 30 m, ...)
    uax::Axis; vax::Axis; sax::Axis
    levels::Vector{Float64}                                # the 29 lower-table levels, km
    mt_rows::Vector{Float64}
    nlon::Int; lstep::Float64
    solar_offset_h::Float64; ho_km::Float64
    K::Vector{Float64}; K_undetermined::Set{Int}           # reference knots (TK x PK, flat), NaN where not fitted
    knee::Float64; slope::Vector{Float64}; slope_undetermined::Set{Int}
end

_floats(x) = x isa AbstractVector && all(v -> v isa Real && !(v isa Bool), x) ? Vector{Float64}(x) : nothing

"""
    WindLayer(payload::AbstractDict, scalars::MarsNearSurfaceScalars.NearSurfaceModel)

Validate the wind part of a deserialized `spaceagora_mars_near_surface_v2` payload against its scalar part, and build
the layer. Stored fields may hold finite values or NaN sentinels, not infinities. The arrays are used as stored; treat
them as read-only.
"""
function WindLayer(d::AbstractDict, s::MarsNearSurfaceScalars.NearSurfaceModel)
    get(d, "format", "") == FORMAT || _payload_error("Wind layer requires payload format $FORMAT.")
    for key in (ARRAY_KEYS..., "wind_meta", "sound_reference", "sound_composition", "wind_rules")
        haskey(d, key) || _payload_error("Near-surface wind payload is missing '$key'.")
    end
    m, ref, comp, rules = d["wind_meta"], d["sound_reference"], d["sound_composition"], d["wind_rules"]
    all(x -> x isa AbstractDict, (m, ref, comp, rules)) || _payload_error("wind_meta, sound_reference, sound_composition and wind_rules must be dictionaries.")
    # The rules this evaluator implements; a payload declaring others is refused.
    get(rules, "clip_fraction", nothing) == CLIP_FRACTION && get(rules, "composition_switch_km", nothing) == SWITCH_KM &&
        get(rules, "switch_band_km", nothing) == SWITCH_BAND_KM ||
        _payload_error("wind_rules must declare clip_fraction $CLIP_FRACTION, composition_switch_km $SWITCH_KM and switch_band_km $SWITCH_BAND_KM.")
    _floats(get(ref, "temperature_levels_K", nothing)) == TK && _floats(get(ref, "pressure_levels_Pa", nothing)) == PK &&
        _floats(get(comp, "slope_levels_K", nothing)) == TK ||
        _payload_error("The sound-speed reference must use the temperature levels $TK and pressure levels $PK.")
    get(comp, "R_u", nothing) == RU || _payload_error("sound_composition R_u must be $RU.")
    levels, lons, mt = _floats(get(m, "levels_km", nothing)), _floats(get(m, "longitudes_deg", nothing)), _floats(get(m, "mtgcm_rows_deg", nothing))
    ulat, vlat, slat = _floats(get(m, "u_knots_deg", nothing)), _floats(get(m, "v_knots_deg", nothing)), _floats(get(m, "s_knots_deg", nothing))
    usides = get(m, "u_sides", nothing)
    any(isnothing, (levels, lons, mt, ulat, vlat, slat)) && _payload_error("wind_meta axes must be real vectors.")
    usides isa AbstractVector{<:Integer} && length(usides) == length(ulat) || _payload_error("wind_meta u_sides must give one integer side per east-wind knot.")
    nl, nlat, nlon = length(levels), s.nlat, s.nlon
    nl >= 2 && levels == s.levels[1:min(nl, end)] && length(s.levels) == nl + 2 ||
        _payload_error("wind_meta levels_km must be the scalar levels without the two MTGCM levels.")
    ho, off = get(m, "ho_km", nothing), get(m, "solar_offset_h", nothing)
    all(x -> x isa Real && !(x isa Bool) && isfinite(x), (ho, off)) || _payload_error("wind_meta ho_km and solar_offset_h must be finite.")
    s.levels[nl+1] == 80.0 + ho && s.levels[nl+2] == 85.0 + ho || _payload_error("The scalar MTGCM levels must be 80 + ho_km and 85 + ho_km.")
    lons == [s.lstep * (j - 1) for j in 1:nlon] || _payload_error("wind_meta longitudes_deg must be the scalar lattice longitudes.")
    slat == [s.lat0 + s.lstep * (i - 1) for i in 1:nlat] || _payload_error("wind_meta s_knots_deg must be the scalar lattice latitudes.")
    for (name, lat) in (("u_knots_deg", ulat), ("v_knots_deg", vlat), ("mtgcm_rows_deg", mt))
        length(lat) >= 2 && all(isfinite, lat) && issorted(lat) || _payload_error("wind_meta $name must be finite and sorted.")
    end
    issorted(mt; lt = <=) && issorted(vlat; lt = <=) || _payload_error("MTGCM rows and north-wind knots must be strictly increasing.")
    sides = Vector{Int}(usides)
    for i in eachindex(ulat)
        sides[i] in (-1, 0, 1) || _payload_error("East-wind knot sides must be -1, 0 or 1.")
        twin = (i > 1 && ulat[i-1] == ulat[i]) || (i < length(ulat) && ulat[i+1] == ulat[i])
        twin == (sides[i] != 0) || _payload_error("East-wind knot $(ulat[i]) must appear twice exactly where it is two-sided.")
        i > 1 && ulat[i-1] == ulat[i] && (sides[i-1], sides[i]) != (-1, 1) && _payload_error("A two-sided knot must list its lower side first.")
    end
    shapes = Dict("wind_level_U" => (nl, length(ulat), nlon), "wind_level_V" => (nl, length(vlat), nlon),
        "wind_mtgcm_U" => (2, length(mt), nlon), "wind_mtgcm_V" => (2, length(mt), nlon),
        "wind_surface_U" => (2, length(ulat), nlon + 1), "wind_surface_V" => (2, length(vlat), nlon),
        "sound_offset_level" => (nl + 2, nlat, nlon), "sound_offset_surface" => (2, nlat, nlon))
    for (key, dims) in shapes
        a = d[key]
        a isa Array{Float64,3} && size(a) == dims || _payload_error("$key must be a Float64 array of size $dims.")
        all(x -> !isinf(x), a) || _payload_error("$key must contain finite values or NaN sentinels.")
    end
    for key in ("wind_level_U", "wind_level_V", "wind_mtgcm_U", "wind_mtgcm_V", "sound_offset_level")
        src = d[key * "_source"]
        src isa Array{UInt8,3} && size(src) == size(d[key]) || _payload_error("$(key)_source must be a UInt8 array of the field's size.")
        all(i -> (src[i] == 0) == isnan(d[key][i]), eachindex(src)) || _payload_error("$key must be NaN exactly where its source code is 0.")
    end
    K = _floats(get(ref, "K", nothing)); slope = _floats(get(comp, "slope_per_kgkmol", nothing)); knee = get(comp, "knee_kgkmol", nothing)
    K !== nothing && length(K) == length(TK) * length(PK) && all(x -> !isinf(x), K) || _payload_error("sound_reference K must hold $(length(TK) * length(PK)) values.")
    slope !== nothing && length(slope) == length(TK) && all(x -> !isinf(x), slope) || _payload_error("sound_composition slope_per_kgkmol must hold $(length(TK)) values.")
    knee isa Real && !(knee isa Bool) && isfinite(knee) || _payload_error("sound_composition knee_kgkmol must be finite.")
    Ku, su = get(ref, "undetermined", nothing), get(comp, "undetermined", nothing)
    Ku isa AbstractVector{<:Integer} && all(in(1:length(K)), Ku) && su isa AbstractVector{<:Integer} && all(in(1:length(TK)), su) ||
        _payload_error("Undetermined knots and slope levels must be valid indices.")
    WindLayer(d["wind_level_U"], d["wind_level_V"], d["wind_mtgcm_U"], d["wind_mtgcm_V"], d["wind_surface_U"], d["wind_surface_V"],
        d["sound_offset_level"], d["sound_offset_surface"], Axis(ulat, sides), Axis(vlat, zeros(Int, length(vlat))),
        Axis(slat, zeros(Int, length(slat))), levels, mt, nlon, s.lstep, Float64(off), Float64(ho),
        K, Set{Int}(Ku), Float64(knee), slope, Set{Int}(su))
end

# ---- horizontal interpolation -----------------------------------------------------------------------------------

"Cell (i, i + 1) and weight of latitude φ on an axis; at a two-sided knot the cell on φ's side, exactly on it the upper side."
function lat_bracket(ax::Axis, φ)
    n = length(ax.lat)
    i = clamp(searchsortedlast(ax.lat, φ), 1, n - 1)
    ax.lat[i+1] > ax.lat[i] || (i = clamp(i + 1, 1, n - 1))
    (i, i + 1, (φ - ax.lat[i]) / (ax.lat[i+1] - ax.lat[i]))
end
"Bracketing stored longitudes (1-based, periodic) and weight; l in [0, 360)."
function lon_bracket(L::WindLayer, l)
    x = mod(l, 360.0) / L.lstep; j = floor(Int, x)
    (j + 1, mod(j + 1, L.nlon) + 1, x - j)
end
"The surface east wind's seam-split longitudes (0+ first, 360- last); exactly at 0 E the west side."
function lon_bracket_seam(L::WindLayer, l)
    l = mod(l, 360.0)
    l == 0.0 && return (L.nlon + 1, L.nlon + 1, 0.0)
    x = l / L.lstep; j = floor(Int, x)
    (j + 1, j + 2, x - j)
end

"Names an unavailable stored node in a refusal."
struct NodeRef
    field::String
    endpoint::String
    lat::Vector{Float64}
    side::Vector{Int}
    seam::Bool
end
function _unavailable(r::NodeRef, i, j, L::WindLayer)
    side = r.side[i] < 0 ? " (lower side)" : r.side[i] > 0 ? " (upper side)" : ""
    lon = r.seam ? (j == 1 ? "0+" : j == L.nlon + 1 ? "360-" : string(L.lstep * (j - 1))) : string(L.lstep * (j - 1))
    throw(DomainError((r.field, r.endpoint, r.lat[i], lon),
        "Near-surface $(r.field) unavailable here: no stored value at $(r.endpoint), latitude knot $(r.lat[i])$(side), longitude $(lon) E"))
end
"Bilinear value of A[k, ·, ·]; nodes with zero weight are skipped, and a needed node that is NaN is refused by name."
function bil(A, k, i1, i2, s, j1, j2, f, r::NodeRef, L::WindLayer)
    v = 0.0
    for (i, wi) in ((i1, 1 - s), (i2, s)), (j, wj) in ((j1, 1 - f), (j2, f))
        w = wi * wj; w == 0 && continue
        a = A[k, i, j]
        isnan(a) && _unavailable(r, i, j, L)
        v += w * a
    end
    v
end

# ---- height rule ---------------------------------------------------------------------------------------------

"An endpoint of the height interpolation: kind :surface (1 = 5 m, 2 = 30 m), :level (lower-table index) or :mtgcm (1 = 80 + ho, 2 = 85 + ho)."
struct Endpoint
    kind::Symbol
    index::Int
end
_name(L::WindLayer, e::Endpoint) = e.kind == :surface ? (e.index == 1 ? "5 m clearance" : "30 m clearance") :
    e.kind == :level ? "the $(L.levels[e.index]) km level" : (e.index == 1 ? "the MTGCM level 80 + ho km" : "the MTGCM level 85 + ho km")

"""
Regime, endpoints and height weight at areoid height z over surface height zs (km). The scalar evaluator has already
refused positions outside the product; at its tolerances (clearance down to 5 m - 1e-9 km) the D5 weight is clamped to 0.
"""
function regime(L::WindLayer, z, zs)
    Δ = z - zs; L1 = max(floor(zs + 0.3) + 1, -5.0)
    Δ <= 0.030 && return (regime = :D5, lo = Endpoint(:surface, 1), hi = Endpoint(:surface, 2), w = max((Δ - 0.005) / 0.025, 0.0))
    if z < L1
        k = findfirst(==(L1), L.levels)
        k === nothing && throw(DomainError(L1, "Near-surface wind unavailable here: first level $(L1) km is not stored"))
        return (regime = :D4, lo = Endpoint(:surface, 2), hi = Endpoint(:level, k), w = (z - zs - 0.03) / (L1 - zs - 0.03))
    end
    if z < 75.0
        i = searchsortedlast(L.levels, z)
        return (regime = :D3, lo = Endpoint(:level, i), hi = Endpoint(:level, i + 1), w = (z - L.levels[i]) / (L.levels[i+1] - L.levels[i]))
    end
    z < 80.0 + L.ho_km && return (regime = :D2, lo = Endpoint(:level, length(L.levels)), hi = Endpoint(:mtgcm, 1), w = (z - 75.0) / (5.0 + L.ho_km))
    (regime = :D1, lo = Endpoint(:mtgcm, 1), hi = Endpoint(:mtgcm, 2), w = (z - 80.0 - L.ho_km) / 5.0)
end

"Pre-clip east or north wind at one endpoint."
function endpoint_wind(L::WindLayer, e::Endpoint, comp::Symbol, φ, l)
    name = comp == :U ? "east wind" : "north wind"
    if e.kind == :mtgcm
        A = comp == :U ? L.MU : L.MV
        k = clamp(searchsortedlast(L.mt_rows, φ), 1, length(L.mt_rows) - 1)
        r = NodeRef(name, _name(L, e), L.mt_rows, zeros(Int, length(L.mt_rows)), false)
        return bil(A, e.index, k, k + 1, (φ - L.mt_rows[k]) / (L.mt_rows[k+1] - L.mt_rows[k]), lon_bracket(L, l)..., r, L)
    end
    ax = comp == :U ? L.uax : L.vax
    i1, i2, s = lat_bracket(ax, φ)
    r = NodeRef(name, _name(L, e), ax.lat, ax.side, e.kind == :surface && comp == :U)
    if e.kind == :surface
        A = comp == :U ? L.SU : L.SV
        return bil(A, e.index, i1, i2, s, (comp == :U ? lon_bracket_seam(L, l) : lon_bracket(L, l))..., r, L)
    end
    bil(comp == :U ? L.U : L.V, e.index, i1, i2, s, lon_bracket(L, l)..., r, L)
end

"The query's stored sound value: the regime's weights; in D2 only the endpoint on the query's side of the switch."
function sound_offset(L::WindLayer, g, φ, l, z)
    i1, i2, s = lat_bracket(L.sax, φ); j1, j2, f = lon_bracket(L, l)
    val(A, k, e) = bil(A, k, i1, i2, s, j1, j2, f, NodeRef("sound speed offset", _name(L, e), L.sax.lat, L.sax.side, false), L)
    endp(e) = e.kind == :surface ? val(L.OS, e.index, e) : e.kind == :level ? val(L.O, e.index, e) : val(L.O, length(L.levels) + e.index, e)
    if g.regime == :D2
        return z <= SWITCH_KM ? val(L.O, length(L.levels), Endpoint(:level, length(L.levels))) :
                                val(L.O, length(L.levels) + 1, Endpoint(:mtgcm, 1))
    end
    w = g.w; (w == 0 ? 0.0 : w * endp(g.hi)) + (w == 1 ? 0.0 : (1 - w) * endp(g.lo))
end

# ---- sound speed ---------------------------------------------------------------------------------------------

pindex(P) = clamp(floor(Int, log10(clamp(P, PK[1], PK[end - 1])) + 2.1), 0, length(PK) - 2)     # native's rule (0-based)
tindex(T) = clamp(floor(Int, clamp(T, TK[1], TK[end - 1]) / 50.0) - 1, 0, length(TK) - 2)

"The reference's knots with nonzero weight at (T, P), in increasing knot order, and their weights."
function cell(T, P)
    i = tindex(T); j = pindex(P)
    t = (clamp(T, TK[1], TK[end]) - TK[i + 1]) / (TK[i + 2] - TK[i + 1])
    p = (clamp(P, PK[1], PK[end]) - PK[j + 1]) / (PK[j + 2] - PK[j + 1])
    k = i * length(PK) + j + 1
    ((k, (1 - t) * (1 - p)), (k + 1, (1 - t) * p), (k + length(PK), t * (1 - p)), (k + length(PK) + 1, t * p))
end
"Reference K at (T, P); refuses a cell the fit did not reach or did not determine."
function Kref(L::WindLayer, T, P)
    total = nothing
    for (k, w) in cell(T, P)
        w == 0 && continue
        (isnan(L.K[k]) || k in L.K_undetermined) &&
            throw(DomainError((T, P), "Near-surface sound speed unavailable here: the reference does not determine the cell of T = $(T) K, P = $(P) Pa"))
        total = total === nothing ? w * L.K[k] : total + w * L.K[k]
    end
    total
end
function Kwater(T)
    t = clamp(T, TW[1], TW[end]); i = clamp(searchsortedlast(TW, t), 1, length(TW) - 1); f = (t - TW[i]) / (TW[i + 1] - TW[i])
    (1 - f) * KW[i] + f * KW[i + 1]
end
"Water mass fraction at native's relative humidity (Savijarvi 1991 saturation form with native's cap); P in Pa."
function qwater(T, P)
    p = P / 100; es = min(p / 1.59, 6.1135 * exp(22.542 * (T - 273.16) / (T + 0.32)))
    RH * 0.407 * es / (p - 0.59 * es)
end
wwater(q) = (xw = MD * q / MW; q * (1 + xw) / (1 + q * (1 + xw)))
Kmix(Kd, w, Kw) = ((1 - w) * Kd + w * Kw * RW) / ((1 - w) + w * RW)
toγ(K) = K / (K - 1)
"Composition term g(M, T): zero at or above the knee; below it (M - knee) times the slope at T."
function composition_term(L::WindLayer, M, T)
    M >= L.knee && return 0.0
    i = tindex(T); t = (clamp(T, TK[1], TK[end]) - TK[i + 1]) / (TK[i + 2] - TK[i + 1])
    total = nothing
    for (k, w) in ((i + 1, 1 - t), (i + 2, t))
        w == 0 && continue
        (isnan(L.slope[k]) || k in L.slope_undetermined) &&
            throw(DomainError((M, T), "Near-surface sound speed unavailable here: the composition slope is not determined at T = $(T) K"))
        total = total === nothing ? w * L.slope[k] : total + w * L.slope[k]
    end
    (M - L.knee) * total
end
"""
γ at a query from its stored value δ: at or below the switch the dry offset with water mixed back in; above it the
remainder plus the composition term at the query's molecular weight RU / R, without water.
"""
function gamma(L::WindLayer, T, P, z, δ, R)
    if z <= SWITCH_KM
        Kd = Kref(L, T, P) + δ
        return toγ(Kmix(Kd, wwater(qwater(T, P)), Kwater(T)))
    end
    toγ(Kref(L, T, P) + composition_term(L, RU / R, T) + δ)
end

# ---- the query -----------------------------------------------------------------------------------------------

"""
    winds(L, S, scalars, lon_east_deg)

Winds and sound speed at a position whose scalar state `scalars` (`MarsNearSurfaceScalars.near_surface_state` of the
scalar model `S`, which also supplies the terrain) is known.
Returns the clipped east, north and vertical winds (m/s), the pre-clip horizontal winds, the sound speed, whether a
component was clipped, the wind regime, the composition side used, and whether the query lies within `SWITCH_BAND_KM`
of the switch.
"""
function winds(L::WindLayer, S::MarsNearSurfaceScalars.NearSurfaceModel, st, lon_east_deg::Real)
    φ, l = st.planetocentric_latitude_deg, mod(Float64(lon_east_deg), 360.0)
    z, zs = st.areoid_height_km, st.surface_height_km
    g = regime(L, z, zs)
    w = g.w
    Ubg = (w == 1 ? 0.0 : (1 - w) * endpoint_wind(L, g.lo, :U, φ, l)) + (w == 0 ? 0.0 : w * endpoint_wind(L, g.hi, :U, φ, l))
    Vbg = (w == 1 ? 0.0 : (1 - w) * endpoint_wind(L, g.lo, :V, φ, l)) + (w == 0 ? 0.0 : w * endpoint_wind(L, g.hi, :V, φ, l))
    t = mod(L.solar_offset_h + l / 15, 24)
    se, sn = slopes(S, φ, l)
    us, vs, W = MarsWindEndpoints.slope_wind(t, φ, z - zs, se, sn, Ubg, Vbg)
    U, V = Ubg + us, Vbg + vs
    δ = sound_offset(L, g, φ, l, z)
    c = sqrt(gamma(L, st.temperature_K, st.pressure_Pa, z, δ, st.gas_constant) * st.pressure_Pa / st.density_kgm3)
    Uc, Vc = clamp(U, -0.7c, 0.7c), clamp(V, -0.7c, 0.7c)
    (wind_east_ms = Uc, wind_north_ms = Vc, wind_up_ms = W, unclipped_wind_east_ms = U, unclipped_wind_north_ms = V,
     sound_speed_ms = c, wind_clipped = Uc != U || Vc != V, wind_regime = g.regime,
     composition_side = z <= SWITCH_KM ? :dry : :switched, in_switch_band = abs(z - SWITCH_KM) <= SWITCH_BAND_KM)
end

"East and north terrain slopes (uncapped), from ±0.25 degree differences of the published terrain."
function slopes(S::MarsNearSurfaceScalars.NearSurfaceModel, φ, l)
    zsf(p, q) = MarsNearSurfaceScalars.terrain(S, S.zs, p, mod(q, 360.0))
    R = 1000MarsNearSurfaceScalars.terrain(S, S.ra, φ, mod(l, 360.0))
    (1000(zsf(φ, l + 0.25) - zsf(φ, l - 0.25)) / (deg2rad(0.5) * cosd(φ) * R), 1000(zsf(φ + 0.25, l) - zsf(φ - 0.25, l)) / (deg2rad(0.5) * R))
end

end # module MarsNearSurfaceWinds
