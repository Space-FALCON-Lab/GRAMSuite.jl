module NearSurfaceAtmosphereTests
using Test, Serialization, SHA
import GRAMSuite
using GRAMSuite: GRAMNearSurfaceAtmosphereModel, near_surface_state, density_state
const NS = GRAMSuite.MarsNearSurfaceScalars

# Synthetic payload with analytic fields (not Mars data): flat terrain at 0.5 km, so L1 = 1 km; linear level states;
# constant surface temperatures; a constant Q_L in every band and 9-degree cell; one volcano block.
const LEVELS = vcat(collect(-5.0:1.0:10.0), collect(15.0:5.0:75.0))
Tlev(z) = 220.0 - 2.0 * z
lnplev(z) = log(700.0) - z / 11.0
const R_GAS, T30, T5, Q = 191.0, 214.0, 216.0, 20.0

function synthetic_payload(; edit! = identity)
    tlat0, tstep, tlon0 = -90.0, 5.0, 0.488
    nr, nc = 37, 72
    zs = fill(0.5, nr, nc); ra = fill(3390.0, nr, nc)
    for i in 1:nr, j in 1:nc          # volcano block around 25 N, 100 E
        lat = tlat0 + tstep * (i - 1); lon = tlon0 + tstep * (j - 1)
        (15.0 <= lat <= 35.0 && 90.0 <= lon <= 110.0) && (zs[i, j] = 12.0)
    end
    lat0, lstep, nlat, nlon = -90.0, 7.5, 25, 48
    nl = length(LEVELS)
    T = [Tlev(LEVELS[a]) for a in 1:nl, i in 1:nlat, j in 1:nlon]
    lnp = [lnplev(LEVELS[a]) for a in 1:nl, i in 1:nlat, j in 1:nlon]
    R = fill(R_GAS, nl, nlat, nlon)
    band = Int[]; cell = Int[]; L = Int[]; status = String[]
    for b in -12:11, c in 0:39
        push!(band, b); push!(cell, c); push!(L, 1); push!(status, "qualified")
    end
    n = length(band)
    d = Dict{String,Any}(
        "format" => "spaceagora_mars_near_surface_scalars_v1", "planet" => "mars", "status" => "synthetic test payload",
        "epoch_utc" => "2001-11-07T11:51:04.794789Z", "radii_km" => (3396.19, 3376.2), "levels_km" => copy(LEVELS),
        "terrain" => Dict{String,Any}("lat0_deg" => tlat0, "lon0_deg" => tlon0, "step_deg" => tstep,
            "surface_height_km" => zs, "areoid_radius_km" => ra),
        "lattice" => Dict{String,Any}("lat0_deg" => lat0, "step_deg" => lstep, "nlat" => nlat, "nlon" => nlon),
        "level_T_K" => T, "level_R" => R, "level_lnp" => lnp, "level_source" => ones(UInt8, nl, nlat, nlon),
        "surface_T30_K" => fill(T30, nlat, nlon), "surface_T5_K" => fill(T5, nlat, nlon),
        "q_models" => Dict{String,Any}("band" => band, "cell" => cell, "L" => L, "order" => fill(1, n),
            "phic_center" => [7.5b + 3.75 for b in band], "lam_center" => [9.0c + 4.5 for c in cell],
            "coef" => hcat(fill(Q, n), zeros(n, 5)), "n_points" => fill(1, n), "status" => status),
        "support" => Dict{String,Any}("zs_refuse_km" => 9.0, "phic_max_deg" => 85.0, "min_clearance_km" => 0.005, "top_areoid_km" => 75.0))
    edit!(d)
    return d
end

function write_payload(d)
    file = tempname() * ".jls"
    open(io -> serialize(io, d), file, "w")
    return file, bytes2hex(open(sha256, file))
end
model_of(d) = (f = write_payload(d); GRAMNearSurfaceAtmosphereModel(surrogate_file = f[1], expected_sha256 = f[2]))

# Geodetic latitude, east longitude and ellipsoid height (m) reaching areoid height z_km at planetocentric latitude phic.
function geodetic_query(m, phic, lam, z_km)
    c = m.core
    r = NS.terrain(c, c.ra, phic, lam) + z_km
    g = phic; h = r - c.a_km
    for _ in 1:50
        rr, pc = NS.radial_position(c.a_km, c.b_km, g, h)
        g += phic - pc; h += r - rr
    end
    return g, lam, 1000h
end
surface_z(m, phic, lam) = NS.terrain(m.core, m.core.zs, phic, lam)
function refusal_message(f)
    try
        f(); return ""
    catch e
        e isa DomainError || rethrow()
        return string(e.msg)
    end
end

const M = model_of(synthetic_payload())

@testset "Near-surface regimes reproduce the analytic rule" begin
    phic, lam = 10.3, 40.7
    zs = surface_z(M, phic, lam); L1 = max(floor(zs + 0.3) + 1, -5.0)
    @test L1 == 1.0
    # D3 between the 3 and 4 km levels
    s = near_surface_state(M, geodetic_query(M, phic, lam, 3.3)...)
    f = (s.areoid_height_km - 3.0)
    @test s.regime == :D3 && s.surface_layer_model_status == "not used" && s.first_level_km == 1.0
    @test s.temperature_K ≈ Tlev(3.0) + (Tlev(4.0) - Tlev(3.0)) * f rtol = 1e-12
    p = exp(lnplev(3.0) + (lnplev(4.0) - lnplev(3.0)) * f)
    @test s.pressure_Pa ≈ p rtol = 1e-12
    @test s.density_kgm3 ≈ p / (R_GAS * s.temperature_K) rtol = 1e-12
    # D4 at 200 m clearance
    s = near_surface_state(M, geodetic_query(M, phic, lam, zs + 0.2)...)
    z, z30 = s.areoid_height_km, s.surface_height_km + 0.03
    @test s.regime == :D4 && s.surface_layer_model_status == "qualified"
    @test s.temperature_K ≈ T30 + (Tlev(1.0) - T30) * (z - z30) / (1.0 - z30) rtol = 1e-12
    H = (Tlev(1.0) + T5) / Q
    @test s.pressure_Pa ≈ exp(lnplev(1.0)) * exp((1.0 - z) / H) rtol = 1e-12
    @test s.gas_constant ≈ R_GAS rtol = 1e-12       # bilinear weights sum to 1 within rounding
    # D5 at 15 m clearance
    s = near_surface_state(M, geodetic_query(M, phic, lam, zs + 0.015)...)
    z = s.areoid_height_km
    @test s.regime == :D5
    @test s.temperature_K ≈ T5 + (T30 - T5) * (z - (s.surface_height_km + 0.005)) / 0.025 rtol = 1e-12
    @test s.pressure_Pa ≈ exp(lnplev(1.0)) * exp((1.0 - z) / H) rtol = 1e-12
end

@testset "Refusals name their reason" begin
    q(phic, lam, z) = () -> near_surface_state(M, geodetic_query(M, phic, lam, z)...)
    @test occursin("Planetocentric latitude beyond", refusal_message(q(86.0, 10.0, 5.0)))
    @test occursin("volcano flanks", refusal_message(q(25.0, 100.0, 20.0)))
    zs = surface_z(M, 10.3, 40.7)
    @test occursin("Below 5 m", refusal_message(q(10.3, 40.7, zs + 0.004)))
    @test occursin("areoid-height top", refusal_message(q(10.3, 40.7, 75.5)))
    @test occursin("Finite geodetic position", refusal_message(() -> near_surface_state(M, NaN, 0.0, 0.0)))
    # an unavailable level at one node refuses the queries that need it, and only those
    mnan = model_of(synthetic_payload(edit! = d -> (d["level_T_K"][7, 13, 7] = NaN)))   # level 1 km at 0 N, 45 E
    @test occursin("component unavailable", refusal_message(() -> near_surface_state(mnan, geodetic_query(mnan, 1.0, 44.0, surface_z(mnan, 1.0, 44.0) + 0.2)...)))
    @test near_surface_state(mnan, geodetic_query(mnan, 1.0, 44.0, 5.5)...).regime == :D3
    msurf = model_of(synthetic_payload(edit! = d -> (d["surface_T30_K"][13, 7] = NaN)))
    @test occursin("surface field unavailable", refusal_message(() -> near_surface_state(msurf, geodetic_query(msurf, 1.0, 44.0, surface_z(msurf, 1.0, 44.0) + 0.2)...)))
    mq = model_of(synthetic_payload(edit! = function (d)
        k = findfirst(i -> d["q_models"]["band"][i] == 2 && d["q_models"]["cell"][i] == 5, eachindex(d["q_models"]["band"]))
        for key in ("band", "cell", "L", "order", "phic_center", "lam_center", "n_points", "status")
            deleteat!(d["q_models"][key], k)
        end
        d["q_models"]["coef"] = d["q_models"]["coef"][setdiff(1:end, k), :]
    end))
    @test occursin("coefficient unavailable", refusal_message(() -> near_surface_state(mq, geodetic_query(mq, 18.0, 49.0, surface_z(mq, 18.0, 49.0) + 0.2)...)))
end

@testset "Longitude just below the terrain origin wraps" begin
    hit = filter(l -> mod(l - M.core.tlon0, 360.0) == 360.0, [prevfloat(0.488, k) for k in 1:1000])
    @test !isempty(hit)
    g, _, h = geodetic_query(M, 20.0, 0.488, 4.0)
    s0 = near_surface_state(M, g, 0.488, h)
    for lam in hit
        @test NS.terrain(M.core, M.core.zs, 20.0, lam) == NS.terrain(M.core, M.core.zs, 20.0, 0.488)
        @test near_surface_state(M, g, lam, h).density_kgm3 ≈ s0.density_kgm3 rtol = 1e-12
    end
end

@testset "Surface-layer model on an edge comes from the neighbour across it" begin
    function edited(which)
        synthetic_payload(edit! = function (d)
            qm = d["q_models"]
            find(b, c) = findfirst(i -> qm["band"][i] == b && qm["cell"][i] == c, eachindex(qm["band"]))
            if which == :band        # own band -1 missing at cell 4; band -2 is a decoy
                own, decoy, across = find(-1, 4), find(-2, 4), find(0, 4)
            else                     # own cell 10 missing in band 1; cell 9 is a decoy
                own, decoy, across = find(1, 10), find(1, 9), find(1, 11)
            end
            qm["status"][decoy] = "decoy two away"; qm["status"][across] = "across the edge"
            for key in ("band", "cell", "L", "order", "phic_center", "lam_center", "n_points", "status")
                deleteat!(qm[key], own)
            end
            qm["coef"] = qm["coef"][setdiff(1:end, own), :]
        end)
    end
    mb = model_of(edited(:band))
    for phic in (-1e-12, 0.0)       # just below and exactly on the equator (band edge 0)
        g, _, h = geodetic_query(mb, phic, 40.0, surface_z(mb, phic, 40.0) + 0.2)
        c = NS.radial_position(mb.core.a_km, mb.core.b_km, g, h / 1000)[2]
        @test abs(c) < 7.5e-9
        @test near_surface_state(mb, g, 40.0, h).surface_layer_model_status == "across the edge"
    end
    mc = model_of(edited(:cell))
    lam = prevfloat(99.0, 4)        # just west of the 99 E cell edge: own cell 10 is missing
    g, _, h = geodetic_query(mc, 10.0, lam, surface_z(mc, 10.0, lam) + 0.2)
    @test near_surface_state(mc, g, lam, h).surface_layer_model_status == "across the edge"
end

@testset "Surface-layer grid corners include the diagonal cell" begin
    function corner_model(keys)
        model_of(synthetic_payload(edit! = function (d)
            fill!(d["terrain"]["surface_height_km"], 0.5)
            qm = d["q_models"]
            keep = findall(i -> (qm["band"][i], qm["cell"][i], qm["L"][i]) in keys, eachindex(qm["band"]))
            qm["status"] = [string((b, c, L)) for (b, c, L) in zip(qm["band"], qm["cell"], qm["L"])]
            for key in ("band", "cell", "L", "order", "phic_center", "lam_center", "n_points", "status")
                qm[key] = qm[key][keep]
            end
            qm["coef"] = qm["coef"][keep, :]
            # Nonconstant Q also checks that the selected cell's coordinates are used.
            qm["coef"][:, 2] .= 0.2
            qm["coef"][:, 4] .= 0.3
        end))
    end
    # Exact corners, both hemispheres, equivalent seam longitudes and all four approach directions.
    points = [(0.0, 0.0), (7.5, 99.0), (-7.5, 81.0), (0.0, 360.0), (0.0, -360.0)]
    append!(points, [(a * 1e-12, b * 1e-12) for a in (-1, 1) for b in (-1, 1)])
    for (phic, lon) in points, clearance in (0.015, 0.2)
        @testset "corner $phic, $lon at $clearance km clearance" begin
            query = geodetic_query(M, phic, lon, 0.5 + clearance)
            actual_phic = NS.radial_position(M.core.a_km, M.core.b_km, query[1], query[3] / 1000)[2]
            # Choose the opposite member of each explicit pair of incident cells.
            eb, ec = round(Int, phic / 7.5), round(Int, mod(lon, 360.0) / 9.0)
            b = actual_phic < 7.5eb ? eb : eb - 1
            c = mod(mod(lon, 360.0) < 9.0ec ? ec : ec - 1, 40)
            key = (b, c, 1)
            m = corner_model([key])
            s = near_surface_state(m, query...)
            @test s.surface_layer_model_status == string(key)
            @test s.regime == (clearance < 0.03 ? :D5 : :D4)
            x = (actual_phic - (7.5b + 3.75)) / 3.75
            y = (mod(lon - (9.0c + 4.5) + 180.0, 360.0) - 180.0) / 4.5
            q = Q + 0.2x + 0.3y
            p = exp(lnplev(1.0)) * exp((1.0 - s.areoid_height_km) * q / (Tlev(1.0) + T5))
            @test s.pressure_Pa ≈ p rtol = 1e-12
            @test s.density_kgm3 ≈ p / (R_GAS * s.temperature_K) rtol = 1e-12
        end
    end
    # Existing primary, latitude-neighbour and longitude-neighbour precedence is preserved.
    keys = [(0, 0, 1), (-1, 0, 1), (0, 39, 1), (-1, 39, 1)]
    for first in eachindex(keys)
        m = corner_model(keys[first:end])
        @test near_surface_state(m, geodetic_query(m, 0.0, 0.0, 0.7)...).surface_layer_model_status == string(keys[first])
    end
    m = corner_model([last(keys)])
    for (phic, lon) in ((1e-6, 0.0), (0.0, 1e-6), (1e-6, 1e-6))
        @test occursin("coefficient unavailable", refusal_message(() -> near_surface_state(m, geodetic_query(m, phic, lon, 0.7)...)))
    end
    empty = corner_model(NTuple{3,Int}[])
    @test occursin("coefficient unavailable", refusal_message(() -> near_surface_state(empty, geodetic_query(empty, 0.0, 0.0, 0.7)...)))
    @test near_surface_state(empty, geodetic_query(empty, 0.0, 0.0, 2.5)...).regime == :D3
end

@testset "Fixed-grid calling convention" begin
    g, _, h = geodetic_query(M, 10.3, 40.7, 2.5)
    s = near_surface_state(M, g, 40.7, h)
    rho, T, wind = density_state(M, h, deg2rad(g), deg2rad(40.7))
    ref = near_surface_state(M, rad2deg(deg2rad(g)), rad2deg(deg2rad(40.7)), h)
    @test rho === ref.density_kgm3 && T === ref.temperature_K
    @test rho ≈ s.density_kgm3 rtol = 1e-12
    @test wind == zeros(3)
    @test density_state(M, h, deg2rad(g), deg2rad(40.7), 1234.5, false) == (rho, T, wind)   # frozen; wind flag has nothing to select
    @test_throws DomainError density_state(M, NaN, 0.1, 0.2)
    @test_throws DomainError density_state(M, h, 0.1, 0.2, Inf)
    @test_throws DomainError density_state(M, 100_000.0, deg2rad(g), deg2rad(40.7))
end

@testset "Construction checks" begin
    d = synthetic_payload()
    file, digest = write_payload(d)
    m = GRAMNearSurfaceAtmosphereModel(surrogate_file = file)
    @test m.source_sha256 == digest && m.core.source_sha256 == digest
    @test !haskey(m.metadata, "level_T_K") && !haskey(m.metadata, "q_models") && m.metadata["support"]["top_areoid_km"] == 75.0
    @test_throws ArgumentError GRAMNearSurfaceAtmosphereModel(surrogate_file = file, expected_sha256 = "0"^64)
    @test_throws ArgumentError GRAMNearSurfaceAtmosphereModel(surrogate_file = file, expected_sha256 = "short")
    @test_throws ArgumentError GRAMNearSurfaceAtmosphereModel(planet = "Venus", surrogate_file = file)
    @test_throws ArgumentError GRAMNearSurfaceAtmosphereModel(surrogate_file = "")
    @test_throws ArgumentError GRAMNearSurfaceAtmosphereModel(surrogate_file = tempname())
    empty = tempname(); touch(empty)
    @test_throws ArgumentError GRAMNearSurfaceAtmosphereModel(surrogate_file = empty)
    lfs = tempname(); write(lfs, "version https://git-lfs.github.com/spec/v1\noid sha256:" * "0"^64 * "\nsize 10\n")
    @test_throws ArgumentError GRAMNearSurfaceAtmosphereModel(surrogate_file = lfs)
    bad(edit!) = GRAMNearSurfaceAtmosphereModel(surrogate_file = write_payload(synthetic_payload(edit! = edit!))[1])
    @test_throws ArgumentError bad(d -> (d["format"] = "spaceagora_gram_static_grid_v1"))
    @test_throws ArgumentError bad(d -> (d["planet"] = "earth"))
    @test_throws ArgumentError bad(d -> delete!(d, "support"))
    @test_throws ArgumentError bad(d -> (d["level_R"] = d["level_R"][:, 1:end-1, :]))
    @test_throws ArgumentError bad(d -> (d["levels_km"] = reverse(d["levels_km"])))
    @test_throws ArgumentError bad(d -> (d["terrain"]["surface_height_km"][1, 1] = NaN))
    @test_throws ArgumentError bad(d -> (d["q_models"]["coef"][1, 1] = Inf))
    @test_throws ArgumentError bad(d -> (d["radii_km"] = (3376.2, 3396.19)))
    notdict = tempname(); open(io -> serialize(io, [1, 2, 3]), notdict, "w")
    @test_throws ArgumentError GRAMNearSurfaceAtmosphereModel(surrogate_file = notdict)
end

# Keep the stored fields aligned while changing the vertical coverage.
function retain_levels!(d, levels)
    indices = [findfirst(==(z), d["levels_km"]) for z in levels]
    d["levels_km"] = d["levels_km"][indices]
    for key in ("level_T_K", "level_R", "level_lnp", "level_source")
        d[key] = d[key][indices, :, :]
    end
    return d
end
flat_payload() = synthetic_payload(edit! = d -> fill!(d["terrain"]["surface_height_km"], 0.5))

@testset "Near-surface payload validation regressions" begin
    @testset "Reject infinities while preserving NaN sentinels" begin
        for key in ("level_T_K", "level_R", "level_lnp", "surface_T30_K", "surface_T5_K")
            for value in (Inf, -Inf)
                d = synthetic_payload()
                d[key][1] = value
                @test_throws ArgumentError model_of(d)
            end
            d = synthetic_payload()
            if ndims(d[key]) == 3
                d[key][findfirst(==(1.0), d["levels_km"]), 13, 7] = NaN
            else
                d[key][13, 7] = NaN
            end
            m = model_of(d)
            @test m isa GRAMNearSurfaceAtmosphereModel
            q = geodetic_query(m, 1.0, 44.0, surface_z(m, 1.0, 44.0) + 0.2)
            @test_throws DomainError near_surface_state(m, q...)
            @test near_surface_state(m, geodetic_query(m, 1.0, 44.0, 5.5)...).regime == :D3
        end
    end

    @testset "Every supported first exposed level is stored exactly" begin
        @test_throws ArgumentError model_of(retain_levels!(flat_payload(), [0.0, 2.0, 5.0, 75.0]))
        # A near match must not be accepted in place of the actual L1.
        d = flat_payload(); d["levels_km"][findfirst(==(1.0), d["levels_km"])] = nextfloat(1.0)
        @test_throws ArgumentError model_of(d)
        # No terrain node selects L1=2, but interpolation between nodes does.
        d = flat_payload(); d["terrain"]["surface_height_km"][:, 2:end] .= 3.5
        retain_levels!(d, [1.0, 3.0, 4.0, 5.0, 75.0])
        @test_throws ArgumentError model_of(d)
        # Sparse higher levels remain valid when the terrain only needs L1=1.
        m = model_of(retain_levels!(flat_payload(), [1.0, 3.0, 5.0, 75.0]))
        @test near_surface_state(m, geodetic_query(m, 10.3, 40.7, 0.7)...).first_level_km == 1.0
        # Unsupported polar heights and volcano heights do not enlarge coverage.
        d = flat_payload(); d["support"]["phic_max_deg"] = 75.0
        d["terrain"]["surface_height_km"][[1, end], :] .= 8.0
        @test model_of(retain_levels!(d, [1.0, 3.0, 75.0])) isa GRAMNearSurfaceAtmosphereModel
        d = synthetic_payload(); d["support"]["zs_refuse_km"] = 8.7
        retain_levels!(d, filter(!=(10.0), d["levels_km"]))
        @test model_of(d) isa GRAMNearSurfaceAtmosphereModel
        # Deep terrain is clamped to the native -5 km first level.
        d = flat_payload(); fill!(d["terrain"]["surface_height_km"], -8.0)
        @test model_of(retain_levels!(d, [-5.0, 1.0, 75.0])) isa GRAMNearSurfaceAtmosphereModel
        @test_throws ArgumentError model_of(retain_levels!(d, [1.0, 75.0]))
    end

    @testset "Advertised top stays within stored levels" begin
        d = synthetic_payload(); d["support"]["top_areoid_km"] = 76.0
        @test_throws ArgumentError model_of(d)
        d["support"]["top_areoid_km"] = 50.0
        m = model_of(d)
        @test near_surface_state(m, geodetic_query(m, 10.3, 40.7, 50.0)...).temperature_K ≈ Tlev(50.0)
        @test_throws DomainError near_surface_state(m, geodetic_query(m, 10.3, 40.7, 50.5)...)
        @test near_surface_state(M, geodetic_query(M, 10.3, 40.7, 75.0)...).temperature_K ≈ Tlev(75.0)
    end
end
end
