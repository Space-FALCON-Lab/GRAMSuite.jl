module NearSurfaceWindsTests
using Test, Serialization, SHA
import GRAMSuite
using GRAMSuite: GRAMNearSurfaceAtmosphereModel, near_surface_state, near_surface_wind_state, near_surface_winds_available, density_state
const NS = GRAMSuite.MarsNearSurfaceScalars
const NW = GRAMSuite.MarsNearSurfaceWinds
const SW = GRAMSuite.MarsWindEndpoints.slope_wind

# Synthetic version 2 payload with analytic fields (not Mars data).
# - Terrain: flat at 0.5 km (L1 = 1 km), except a ridge between 195.488 and 295.488 E. It rises eastward by 0.01 km per
#   degree up to 245.488 E, then falls back.
# - Scalars: T linear in height, ln p linear, R = 192 J/(kg K), so M = R_u / R = 43.30 kg/kmol (below the 43.5 knee).
# - Winds: linear in latitude and uniform in longitude, with jumps of +3 m/s at 78.75 N and +2 m/s more at 82.5 N in
#   the east wind.
# - The surface east wind's 360- column is 1 m/s above its 0+ column.
# - Sound speed: the reference K is 3.5 everywhere, and offsets are 0.001 times the level index.
const HO = 0.0323
const LOW = vcat(collect(-5.0:1.0:10.0), collect(15.0:5.0:75.0))
const LEVELS = vcat(LOW, [80.0 + HO, 85.0 + HO])
Tlev(z) = 240.0 - z
lnplev(z) = log(700.0) - z / 11.0
const R_GAS, T30, T5, Q, OFF = 192.0, 214.0, 216.0, 20.0, 12.0
jump(φ, side) = (φ > 78.75 || (φ == 78.75 && side > 0) ? 3.0 : 0.0) + (φ > 82.5 || (φ == 82.5 && side > 0) ? 2.0 : 0.0)
jumpq(φ) = (φ >= 78.75 ? 3.0 : 0.0) + (φ >= 82.5 ? 2.0 : 0.0)           # a query exactly on a jump takes the upper side
const ULAT, USIDE = let ks = collect(-90.0:5.0:90.0), sd = zeros(Int, 37)
    append!(ks, [78.75, 78.75, 82.5, 82.5]); append!(sd, [-1, 1, -1, 1])
    o = sortperm(collect(zip(ks, sd))); ks[o], sd[o]
end
const VLAT, MT = collect(-90.0:5.0:90.0), collect(-87.5:5.0:87.5)
fU(k, φ) = 10.0 + 0.5k + 0.1φ                  # lower-table east wind at level index k, before jumps
fV(k, φ) = -5.0 + 0.2k - 0.05φ
fMU(m, φ) = 30.0 + 5m + 0.1φ; fMV(m, φ) = -10.0 + m - 0.05φ
fSU(m, φ) = 2.0 + m + 0.1φ; fSV(m, φ) = 1.0 + m - 0.05φ
const SLOPES = [NaN, 0.023, 0.028, 0.041, NaN, NaN, NaN]

function payload(; format = "spaceagora_mars_near_surface_v2", winds = true, lattice_step = 7.5, longitude_count = 48, edit! = identity)
    tlat0, tstep, tlon0 = -90.0, 5.0, 0.488
    nr, nc = 37, 72
    zs = [0.5 + 0.01 * max(0.0, 50.0 - abs(tlon0 + tstep * (j - 1) - 245.488)) for i in 1:nr, j in 1:nc]
    ra = fill(3390.0, nr, nc)
    lat0, lstep = -90.0, lattice_step
    nlat, nlon = ceil(Int, 180 / lstep) + 1, longitude_count
    nl = length(LEVELS)
    band = Int[]; cell = Int[]; L = Int[]
    for b in -12:11, c in 0:39, l in (1, 2)
        push!(band, b); push!(cell, c); push!(L, l)
    end
    n = length(band)
    d = Dict{String,Any}(
        "format" => format, "planet" => "mars", "status" => "synthetic test payload",
        "epoch_utc" => "2001-11-07T11:51:04.794789Z", "radii_km" => (3396.19, 3376.2), "levels_km" => copy(winds ? LEVELS : LOW),
        "terrain" => Dict{String,Any}("lat0_deg" => tlat0, "lon0_deg" => tlon0, "step_deg" => tstep, "surface_height_km" => zs, "areoid_radius_km" => ra),
        "lattice" => Dict{String,Any}("lat0_deg" => lat0, "step_deg" => lstep, "nlat" => nlat, "nlon" => nlon),
        "level_T_K" => [Tlev(z) for z in (winds ? LEVELS : LOW), i in 1:nlat, j in 1:nlon],
        "level_R" => fill(R_GAS, winds ? nl : 29, nlat, nlon),
        "level_lnp" => [lnplev(z) for z in (winds ? LEVELS : LOW), i in 1:nlat, j in 1:nlon],
        "level_source" => ones(UInt8, winds ? nl : 29, nlat, nlon),
        "surface_T30_K" => fill(T30, nlat, nlon), "surface_T5_K" => fill(T5, nlat, nlon),
        "q_models" => Dict{String,Any}("band" => band, "cell" => cell, "L" => L, "order" => fill(1, n),
            "phic_center" => [7.5b + 3.75 for b in band], "lam_center" => [9.0c + 4.5 for c in cell],
            "coef" => hcat(fill(Q, n), zeros(n, 5)), "n_points" => fill(1, n), "status" => fill("qualified", n)),
        "support" => Dict{String,Any}("zs_refuse_km" => 9.0, "phic_max_deg" => 85.0, "min_clearance_km" => 0.005, "top_areoid_km" => winds ? 81.0 : 75.0))
    if winds
        nu, nv = length(ULAT), length(VLAT)
        d["wind_level_U"] = [fU(k, ULAT[i]) + jump(ULAT[i], USIDE[i]) for k in 1:29, i in 1:nu, j in 1:nlon]
        d["wind_level_V"] = [fV(k, VLAT[i]) for k in 1:29, i in 1:nv, j in 1:nlon]
        d["wind_mtgcm_U"] = [fMU(m, MT[i]) for m in 1:2, i in 1:36, j in 1:nlon]
        d["wind_mtgcm_V"] = [fMV(m, MT[i]) for m in 1:2, i in 1:36, j in 1:nlon]
        d["wind_surface_U"] = [fSU(m, ULAT[i]) + jump(ULAT[i], USIDE[i]) + (j == nlon + 1 ? 1.0 : 0.0) for m in 1:2, i in 1:nu, j in 1:nlon+1]
        d["wind_surface_V"] = [fSV(m, VLAT[i]) for m in 1:2, i in 1:nv, j in 1:nlon]
        d["sound_offset_level"] = [0.001k for k in 1:nl, i in 1:nlat, j in 1:nlon]
        d["sound_offset_surface"] = [m == 1 ? 0.0005 : 0.0007 for m in 1:2, i in 1:nlat, j in 1:nlon]
        for key in ("wind_level_U", "wind_level_V", "wind_mtgcm_U", "wind_mtgcm_V", "sound_offset_level")
            d[key * "_source"] = ones(UInt8, size(d[key]))
        end
        d["wind_meta"] = Dict{String,Any}("levels_km" => copy(LOW), "longitudes_deg" => [lstep * (j - 1) for j in 1:nlon],
            "s_knots_deg" => [lat0 + lstep * (i - 1) for i in 1:nlat], "u_knots_deg" => copy(ULAT), "u_sides" => copy(USIDE),
            "v_knots_deg" => copy(VLAT), "mtgcm_rows_deg" => copy(MT), "ho_km" => HO, "solar_offset_h" => OFF)
        d["sound_reference"] = Dict{String,Any}("K" => fill(3.5, 56), "undetermined" => Int[], "nrows" => 1,
            "temperature_levels_K" => copy(NW.TK), "pressure_levels_Pa" => copy(NW.PK))
        d["sound_composition"] = Dict{String,Any}("knee_kgkmol" => 43.5, "slope_per_kgkmol" => copy(SLOPES), "slope_levels_K" => copy(NW.TK),
            "undetermined" => Int[], "R_u" => NW.RU)
        d["wind_rules"] = Dict{String,Any}("clip_fraction" => 0.7, "composition_switch_km" => 80.0, "switch_band_km" => 1e-9)
    end
    edit!(d)
    d
end
function model_of(d)
    file = tempname() * ".jls"; open(io -> serialize(io, d), file, "w")
    GRAMNearSurfaceAtmosphereModel(surrogate_file = file, expected_sha256 = bytes2hex(open(sha256, file)))
end
"Geodetic latitude, east longitude and ellipsoid height (m) reaching areoid height z (km) at planetocentric latitude phic."
function at(m, phic, lam, z)
    c = m.core; r = NS.terrain(c, c.ra, phic, lam) + z; g = phic; h = r - c.a_km
    for _ in 1:60
        rr, pc = NS.radial_position(c.a_km, c.b_km, g, h); g += phic - pc; h += r - rr
    end
    g, lam, 1000h
end
zsurf(m, phic, lam) = NS.terrain(m.core, m.core.zs, phic, lam)
query(m, phic, lam, z) = near_surface_wind_state(m, at(m, phic, lam, z)...)
message(f) = try; f(); ""; catch e; e isa DomainError ? string(e.msg) : rethrow(); end

const M = model_of(payload())
const M1 = model_of(payload(format = "spaceagora_mars_near_surface_scalars_v1", winds = false))

@testset "Wind longitude lattice covers one full period" begin
    # Matching scalar/wind axes and array shapes are not sufficient: the stored
    # columns must span one 360-degree period, without a gap or duplicate wrap.
    for step in (5.0, 8.0, 7.500001)
        @test_throws ArgumentError model_of(payload(lattice_step = step))
    end
    for step in (prevfloat(7.5), 7.5, nextfloat(7.5))
        roundoff = model_of(payload(lattice_step = step))
        L = roundoff.winds
        i, j, w = NW.lon_bracket(L, prevfloat(360.0))
        @test (i, j) == (L.nlon, 1)
        @test 0.0 <= w <= 1.0
        i, j, w = NW.lon_bracket_seam(L, prevfloat(360.0))
        @test (i, j) == (L.nlon, L.nlon + 1)
        @test 0.0 <= w <= 1.0
        @test NW.lon_bracket_seam(L, 360.0) == (L.nlon + 1, L.nlon + 1, 0.0)
        for z in (0.52, 2.0, 80.5)
            @test isfinite(query(roundoff, 10.3, prevfloat(360.0), z).wind_east_ms)
        end
    end
    fine = model_of(payload(lattice_step = 5.0, longitude_count = 72))
    for m in (M, fine), z in (0.52, 2.0, 80.5), longitude in (0.0, 300.0, 359.999)
        w = query(m, 10.3, longitude, z)
        @test all(isfinite, (w.wind_east_ms, w.wind_north_ms, w.wind_up_ms, w.sound_speed_ms))
        for wrapped in (longitude - 360, longitude + 360)
            other = query(m, 10.3, wrapped, z)
            @test other.wind_east_ms ≈ w.wind_east_ms
            @test other.wind_north_ms ≈ w.wind_north_ms
            @test other.sound_speed_ms ≈ w.sound_speed_ms
        end
    end
end

@testset "Version 1 payloads keep zero wind" begin
    @test !near_surface_winds_available(M1) && near_surface_winds_available(M)
    g, l, h = at(M1, 10.3, 40.7, 2.0)
    @test density_state(M1, h, deg2rad(g), deg2rad(l))[3] == zeros(3)
    @test_throws ArgumentError near_surface_wind_state(M1, g, l, h)
    @test M1.winds === nothing && !haskey(M1.metadata, "wind_meta")
    @test !any(k -> haskey(M.metadata, k), GRAMSuite.MarsNearSurfaceScalars.WIND_ARRAY_KEYS)   # stored fields are not metadata
end

@testset "Scalars are the scalar evaluator's" begin
    for (phic, lam, z) in ((10.3, 40.7, 2.0), (-40.0, 300.0, 0.52), (60.0, 120.0, 50.0), (22.5, 200.0, 80.5))
        g, l, h = at(M, phic, lam, z)
        s, w = near_surface_state(M, g, l, h), near_surface_wind_state(M, g, l, h)
        @test all(k -> getfield(w, k) === getfield(s, k), keys(s))
    end
end

@testset "Regimes and height weights" begin
    phic, lam = 10.3, 40.7
    zs = zsurf(M, phic, lam)
    # D5: the 5 m and 30 m surface fields
    w = query(M, phic, lam, zs + 0.0175); φ = w.planetocentric_latitude_deg; a = (w.areoid_height_km - w.surface_height_km - 0.005) / 0.025
    @test w.wind_regime == :D5
    @test w.unclipped_wind_east_ms ≈ (1 - a) * fSU(1, φ) + a * fSU(2, φ) rtol = 1e-12
    @test w.unclipped_wind_north_ms ≈ (1 - a) * fSV(1, φ) + a * fSV(2, φ) rtol = 1e-12
    # D4: the 30 m field and the first level (1 km, index 7)
    w = query(M, phic, lam, zs + 0.3); φ = w.planetocentric_latitude_deg; z = w.areoid_height_km
    a = (z - w.surface_height_km - 0.03) / (1.0 - w.surface_height_km - 0.03)
    @test w.wind_regime == :D4
    @test w.unclipped_wind_east_ms ≈ (1 - a) * fSU(2, φ) + a * fU(7, φ) rtol = 1e-12
    # D3: the 10 and 15 km levels (indices 16 and 17)
    w = query(M, phic, lam, 12.5); φ = w.planetocentric_latitude_deg; a = (w.areoid_height_km - 10.0) / 5.0
    @test w.wind_regime == :D3
    @test w.unclipped_wind_east_ms ≈ (1 - a) * fU(16, φ) + a * fU(17, φ) rtol = 1e-12
    @test w.unclipped_wind_north_ms ≈ (1 - a) * fV(16, φ) + a * fV(17, φ) rtol = 1e-12
    # D2: the 75 km level and 80 + ho
    w = query(M, phic, lam, 78.0); φ = w.planetocentric_latitude_deg; a = (w.areoid_height_km - 75.0) / (5.0 + HO)
    @test w.wind_regime == :D2
    @test w.unclipped_wind_east_ms ≈ (1 - a) * fU(29, φ) + a * fMU(1, φ) rtol = 1e-12
    # D1: the two MTGCM levels
    w = query(M, phic, lam, 80.8); φ = w.planetocentric_latitude_deg; a = (w.areoid_height_km - 80.0 - HO) / 5.0
    @test w.wind_regime == :D1
    @test w.unclipped_wind_east_ms ≈ (1 - a) * fMU(1, φ) + a * fMU(2, φ) rtol = 1e-12
    @test w.unclipped_wind_north_ms ≈ (1 - a) * fMV(1, φ) + a * fMV(2, φ) rtol = 1e-12
    # above the slope layer no vertical wind; on flat terrain only roundoff of the interpolated surface
    @test w.wind_up_ms == 0.0 && abs(query(M, phic, lam, zs + 0.3).wind_up_ms) < 1e-12
    # the 1.1.0 tolerance zone: 1.1.0 serves clearance down to 5 m - 1e-9 km, and the D5 weight is clamped at 0
    w = query(M, phic, lam, zs + 0.005 - 4e-10)
    @test w.wind_regime == :D5
    @test w.unclipped_wind_east_ms ≈ fSU(1, w.planetocentric_latitude_deg) rtol = 1e-12
end

@testset "Two-sided east-wind knots and the 0 E seam" begin
    for (phic, expect) in ((78.45, 0.0), (79.05, 3.0), (82.2, 3.0), (82.8, 5.0))
        w = query(M, phic, 40.7, 12.5); φ = w.planetocentric_latitude_deg; a = (w.areoid_height_km - 10.0) / 5.0
        @test jumpq(φ) == expect
        @test w.unclipped_wind_east_ms ≈ (1 - a) * fU(16, φ) + a * fU(17, φ) + expect rtol = 1e-12
        @test w.unclipped_wind_north_ms ≈ (1 - a) * fV(16, φ) + a * fV(17, φ) rtol = 1e-12   # the north wind has no jump
    end
    # exactly at 0 E the surface east wind takes the 360- column (west side); just east of it, the 0+ column
    zs = zsurf(M, 10.3, 0.0)
    g, _, h = at(M, 10.3, 0.0, zs + 0.0175)
    w0 = near_surface_wind_state(M, g, 0.0, h); φ = w0.planetocentric_latitude_deg; a = (w0.areoid_height_km - w0.surface_height_km - 0.005) / 0.025
    @test w0.unclipped_wind_east_ms ≈ (1 - a) * fSU(1, φ) + a * fSU(2, φ) + 1.0 rtol = 1e-12
    w1 = near_surface_wind_state(M, g, 1e-6, h)
    @test w1.unclipped_wind_east_ms ≈ (1 - a) * fSU(1, w1.planetocentric_latitude_deg) + a * fSU(2, w1.planetocentric_latitude_deg) rtol = 1e-9
end

@testset "Slope and vertical winds on tilted terrain" begin
    phic, lam = 10.0, 220.0
    zs = zsurf(M, phic, lam); @test zs ≈ 0.5 + 0.01 * (220.0 - 195.488) rtol = 1e-12
    w = query(M, phic, lam, zs + 0.5)
    φ, z, zq = w.planetocentric_latitude_deg, w.areoid_height_km, w.surface_height_km
    a = (z - zq - 0.03) / (2.0 - zq - 0.03)                               # D4 up to L1 = 2 km (index 8)
    Ubg = (1 - a) * fSU(2, φ) + a * fU(8, φ); Vbg = (1 - a) * fSV(2, φ) + a * fV(8, φ)
    se = 1000 * 0.01 / (deg2rad(1.0) * cosd(φ) * 1000 * 3390.0)          # the terrain's east slope, exact
    t = mod(OFF + lam / 15, 24)
    us, vs, W = SW(t, φ, z - zq, se, 0.0, Ubg, Vbg)
    @test w.wind_regime == :D4 && W != 0.0 && us != 0.0
    @test w.unclipped_wind_east_ms ≈ Ubg + us rtol = 1e-9
    @test w.unclipped_wind_north_ms ≈ Vbg + vs rtol = 1e-9
    @test w.wind_up_ms ≈ W rtol = 1e-9
end

@testset "Sound speed, the composition switch and clipping" begin
    phic, lam = 10.3, 40.7
    # dry side with water at or below 80.0 km (D2 uses the 75 km offset there, index 29)
    w = query(M, phic, lam, 79.9)
    γ = NW.toγ(NW.Kmix(3.5 + 0.029, NW.wwater(NW.qwater(w.temperature_K, w.pressure_Pa)), NW.Kwater(w.temperature_K)))
    @test w.composition_side == :dry && !w.in_switch_band
    @test w.sound_speed_ms ≈ sqrt(γ * w.pressure_Pa / w.density_kgm3) rtol = 1e-12
    @test γ != 3.529 / 2.529                                              # water changes the dry value
    # switched side above 80.0 km in D2: the 80 + ho remainder (index 30) plus the composition at M = R_u / R
    for z in (80.02, 80.8)
        w = query(M, phic, lam, z); T = w.temperature_K
        t = (T - 150.0) / 50.0; g = (NW.RU / R_GAS - 43.5) * ((1 - t) * SLOPES[3] + t * SLOPES[4])
        δ = z < 80.0 + HO ? 0.030 : (1 - (w.areoid_height_km - 80.0 - HO) / 5.0) * 0.030 + (w.areoid_height_km - 80.0 - HO) / 5.0 * 0.031
        K = 3.5 + g + δ
        @test w.composition_side == :switched && 150 < T < 200
        @test w.sound_speed_ms ≈ sqrt(K / (K - 1) * w.pressure_Pa / w.density_kgm3) rtol = 1e-12
    end
    # the side follows the query's own height; within 1e-9 km of 80.0 km it may differ from native's, and is flagged
    w = query(M, phic, lam, 80.0)
    @test w.in_switch_band && abs(w.areoid_height_km - 80.0) <= 1e-9
    @test w.composition_side == (w.areoid_height_km <= 80.0 ? :dry : :switched)
    # band membership by the evaluator's own height, both edges pinned to the last double inside
    for sgn in (1, -1)
        z = 80.0 + sgn * 1e-9
        while !NW.in_switch_band(z); z = sgn > 0 ? prevfloat(z) : nextfloat(z); end
        while NW.in_switch_band(sgn > 0 ? nextfloat(z) : prevfloat(z)); z = sgn > 0 ? nextfloat(z) : prevfloat(z); end
        @test NW.in_switch_band(z) && !NW.in_switch_band(sgn > 0 ? nextfloat(z) : prevfloat(z))
        @test abs(abs(z - 80.0) - 1e-9) < 1e-13                             # the edge is at 1e-9 km, to roundoff
    end
    @test NW.in_switch_band(80.0) && !NW.in_switch_band(80.0 + 2e-9) && !NW.in_switch_band(80.0 - 2e-9)
    # clipping: horizontal components at 0.7c, the vertical wind never
    Mc = model_of(payload(edit! = d -> (d["wind_level_U"][16:17, :, :] .= 400.0)))
    w = query(Mc, phic, lam, 12.5)
    @test w.unclipped_wind_east_ms == 400.0 && w.wind_east_ms == 0.7 * w.sound_speed_ms && w.wind_clipped
    @test w.east_clipped && !w.north_clipped && w.wind_north_ms == w.unclipped_wind_north_ms
    Mneg = model_of(payload(edit! = d -> (d["wind_level_V"][16:17, :, :] .= -400.0)))
    wn = query(Mneg, phic, lam, 12.5)
    @test wn.north_clipped && !wn.east_clipped && wn.wind_north_ms == -0.7 * wn.sound_speed_ms && wn.unclipped_wind_north_ms == -400.0
    @test !query(M, phic, lam, 12.5).wind_clipped
    # the limit exactly: east and north, both signs, just below, at and just above 0.7c. The sound speed does not depend
    # on the stored winds, so the limit is the reference query's; the stored levels are set so that the unclipped wind
    # at the query is exactly the target.
    limit = 0.7 * query(M, phic, lam, 12.5).sound_speed_ms
    function exactly(component, target)
        key, field = component === :east ? ("wind_level_U", :unclipped_wind_east_ms) : ("wind_level_V", :unclipped_wind_north_ms)
        stored = target
        for _ in 1:64
            m = model_of(payload(edit! = d -> (d[key][16:17, :, :] .= stored)))
            got = getproperty(query(m, phic, lam, 12.5), field)
            got == target && return m
            stored = got < target ? nextfloat(stored) : prevfloat(stored)
        end
        error("No stored level gives exactly $target")
    end
    for component in (:east, :north), sgn in (1.0, -1.0), (case, target) in ((:below, prevfloat(limit)), (:at, limit), (:above, nextfloat(limit)))
        q = query(exactly(component, sgn * target), phic, lam, 12.5)
        u, out, flag, other = component === :east ? (q.unclipped_wind_east_ms, q.wind_east_ms, q.east_clipped, q.north_clipped) :
                                                    (q.unclipped_wind_north_ms, q.wind_north_ms, q.north_clipped, q.east_clipped)
        @test 0.7 * q.sound_speed_ms == limit && u == sgn * target
        @test flag == (case !== :below) && !other && q.wind_clipped == flag
        @test out == (case === :below ? u : sgn * limit)
    end
    # a component counts as clipped when its unclipped wind reaches the limit (as in W2's classification)
    for (m, z) in ((M, 12.5), (Mc, 12.5), (Mneg, 12.5), (M, 80.8), (M, zsurf(M, phic, lam) + 0.0175))
        q = query(m, phic, lam, z); limit = 0.7 * q.sound_speed_ms
        @test q.east_clipped == (abs(q.unclipped_wind_east_ms) >= limit) && q.north_clipped == (abs(q.unclipped_wind_north_ms) >= limit)
        @test q.wind_clipped == (q.east_clipped || q.north_clipped) && q.wind_east_ms == clamp(q.unclipped_wind_east_ms, -limit, limit)
    end
    # the stored sound value used: in D3 the levels' offsets (0.001 times the level index) by the height weight
    q = query(M, phic, lam, 12.5); a = (q.areoid_height_km - 10.0) / 5.0
    @test q.sound_offset ≈ (1 - a) * 0.016 + a * 0.017 rtol = 1e-12
    # the fixed-grid convention returns the clipped winds
    g, l, h = at(Mc, phic, lam, 12.5)
    ref = near_surface_wind_state(Mc, rad2deg(deg2rad(g)), rad2deg(deg2rad(l)), h)
    @test density_state(Mc, h, deg2rad(g), deg2rad(l))[3] == [ref.wind_east_ms, ref.wind_north_ms, ref.wind_up_ms]
    @test density_state(Mc, h, deg2rad(g), deg2rad(l), 0.0, false) == density_state(Mc, h, deg2rad(g), deg2rad(l))
end

@testset "Refusals name what is unavailable" begin
    phic, lam = 10.3, 40.7
    # an unavailable east-wind node at the 15 km level, knot 10.0, longitude 37.5 (index 6)
    i = findfirst(==(10.0), ULAT)
    Mn = model_of(payload(edit! = d -> (d["wind_level_U"][17, i, 6] = NaN; d["wind_level_U_source"][17, i, 6] = 0x00)))
    msg = message(() -> query(Mn, phic, lam, 12.5))
    @test occursin("east wind", msg) && occursin("15.0 km level", msg) && occursin("latitude knot 10.0", msg) && occursin("longitude 37.5", msg)
    @test query(Mn, phic, lam, 5.5).wind_regime == :D3                    # a query that does not need it is served
    # a reference cell the fit did not determine: T in (200, 250) K, P in (100, 1000) Pa, knot 3 * 8 + 4 + 1
    Mk = model_of(payload(edit! = d -> (d["sound_reference"]["K"][29] = NaN)))
    @test occursin("reference does not determine", message(() -> query(Mk, phic, lam, 12.5)))
    @test query(Mk, phic, lam, 50.0).sound_speed_ms > 0                     # another cell is served
    # a composition slope not determined at the query's temperature (150 K level)
    Ms = model_of(payload(edit! = d -> (d["sound_composition"]["slope_per_kgkmol"][3] = NaN)))
    @test occursin("composition slope is not determined", message(() -> query(Ms, phic, lam, 80.8)))
    @test query(Ms, phic, lam, 12.5).composition_side == :dry               # no composition at or below the switch
    # scalar refusals come first, unchanged
    @test occursin("Below 5 m", message(() -> query(M, phic, lam, zsurf(M, phic, lam) + 0.004)))
end

@testset "Payload validation" begin
    bad(edit!) = @test_throws ArgumentError model_of(payload(edit! = edit!))
    bad(d -> delete!(d, "wind_meta"))
    bad(d -> delete!(d, "wind_rules"))
    bad(d -> (d["wind_rules"]["clip_fraction"] = 0.8))
    bad(d -> (d["wind_level_V"] = d["wind_level_V"][:, 1:end-1, :]))
    bad(d -> (d["wind_level_U"][3, 4, 5] = NaN))                           # NaN where the source code is not 0
    bad(d -> (d["wind_surface_U"][1, 2, 3] = Inf))
    bad(d -> (d["wind_meta"]["u_sides"][findfirst(==(1), d["wind_meta"]["u_sides"])] = 0))
    bad(d -> (d["wind_meta"]["ho_km"] = 0.03))                             # the MTGCM levels must be 80 + ho and 85 + ho
    bad(d -> (d["sound_reference"]["temperature_levels_K"] = collect(50.0:50.0:350.0) .+ 1))
    bad(d -> (d["sound_composition"]["R_u"] = 8314.0))
    bad(d -> (d["levels_km"][29] = 74.0; d["wind_meta"]["levels_km"][29] = 74.0))   # the rule's D3 ends at 75 km
    bad(d -> (d["support"]["min_clearance_km"] = 0.01))                   # the wind rule's 5 m endpoint
    bad(d -> (d["support"]["top_areoid_km"] = 80.0))                       # top below the first MTGCM level
    bad(d -> (d["wind_meta"]["mtgcm_rows_deg"] = collect(-82.5:5.0:82.5)))   # an axis short of the supported latitudes
    # version 1 semantics are unchanged: under the version 1 format the extra keys are not a wind layer
    @test model_of(payload(edit! = d -> (d["format"] = "spaceagora_mars_near_surface_scalars_v1"))).winds === nothing
end
end
