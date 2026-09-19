module GridAtmosphereTests
using Test
using Serialization
using SHA

import GRAMSuite
const GridAPI = GRAMSuite

# This polynomial is linear in each coordinate separately, including cross
# terms. Trilinear interpolation must reproduce it on nonuniform cells.
# The coordinates here are metres and degrees, independently of the reader's
# conversion to radians. Every synthetic value is positive where required.
function analytic_fields(h_m, latitude_deg, longitude_deg)
    h, a, b = Float64(h_m), Float64(latitude_deg), Float64(longitude_deg)
    return (
        1.0 + 0.002h + 0.003a + 0.001b + 1e-7h*a + 1e-6a*b + 1e-8h*a*b,
        300.0 + 0.01h + 0.1a + 0.05b + 1e-5h*a - 1e-5a*b,
        -10.0 + 0.03h + 0.2a - 0.01b + 1e-6h*a*b,
        4.0 - 0.02h + 0.05a + 0.02b - 1e-5h*b,
        0.1 + 0.001h - 0.01a + 0.003b + 1e-6a*b,
    )
end

const SYNTHETIC_FIELD_KEYS = (
    "density_kgm3", "temperature_K", "wind_ew_ms", "wind_ns_ms", "wind_up_ms",
)

function synthetic_payload(; planet="mars")
    altitude_km = [0.0, 0.25, 1.0]
    latitude_deg = [-60.0, -10.0, 70.0]
    longitude_deg = [0.0, 80.0, 210.0, 300.0]
    fields = Dict{String,Any}()
    for (component, key) in enumerate(SYNTHETIC_FIELD_KEYS)
        fields[key] = [
            analytic_fields(h * 1000, a, b)[component]
            for h in altitude_km, a in latitude_deg, b in longitude_deg
        ]
    end
    return Dict{String,Any}(
        "status" => "ok",
        "type" => "surrogate_trilinear",
        "planet" => planet,
        "grid" => Dict{String,Any}(
            "alt_km" => altitude_km,
            "lat_deg" => latitude_deg,
            "lon_deg" => longitude_deg,
        ),
        "fields" => fields,
        "wind_mode" => "nominal",
        "monte_carlo" => "off",
        "dust" => "synthetic",
        "provenance" => "Analytic synthetic test fixture, not atmospheric data",
    )
end

function write_payload(path, payload=synthetic_payload())
    open(path, "w") do io
        serialize(io, payload)
    end
    return path
end

function check_state(actual, expected; rtol=2e-13, atol=2e-13)
    @test length(actual) == 3
    flattened = (actual[1], actual[2], Tuple(actual[3])...)
    @test length(flattened) == 5
    for (observed, desired) in zip(flattened, expected)
        @test isapprox(observed, desired; rtol=rtol, atol=atol)
    end
end

function argument_error_text(f)
    caught = try
        f()
        nothing
    catch err
        err
    end
    @test caught isa ArgumentError
    return caught isa Exception ? sprint(showerror, caught) : ""
end

function construct_from_payload(directory, payload)
    path = write_payload(joinpath(directory, "case_surrogate.jls"), payload)
    return GridAPI.GRAMGridAtmosphereModel(; planet="mars", surrogate_file=path)
end

@testset "Native-free grid atmosphere" begin
    mktempdir() do directory
        source = write_payload(joinpath(directory, "mars_surrogate.jls"))
        actual_sha = bytes2hex(sha256(read(source)))
        native_wrapper_before = GridAPI._GRAM_WRAPPER[]
        native_wrapper_file_before = GridAPI._GRAM_WRAPPER_FILE[]
        @test native_wrapper_before === nothing

        withenv("GRAM_ROOT" => joinpath(directory, "absent_native_root"),
                "GRAM_LIB" => joinpath(directory, "absent_native_library")) do
            model = GridAPI.GRAMGridAtmosphereModel(;
                planet="mars", surrogate_file=source, expected_sha256=actual_sha,
            )

            @testset "Schema, units, metadata and source identity" begin
                @test model.surrogate isa GridAPI.GRAMOfflineSurrogate
                @test model.surrogate.planet_name == "mars"
                @test model.surrogate.source_file == normpath(source)
                @test model.source_sha256 == actual_sha
                @test model.above_grid === :error
                @test model.vacuum_temperature == 200.0
                @test model.surrogate.alt_nodes_m == [0.0, 250.0, 1000.0]
                @test model.surrogate.lat_nodes_rad ≈ deg2rad.([-60.0, -10.0, 70.0])
                @test model.surrogate.lon_nodes_rad ≈ deg2rad.([0.0, 80.0, 210.0, 300.0])
                expected_metadata = Dict(k => v for (k, v) in synthetic_payload()
                                         if k != "grid" && k != "fields")
                @test model.metadata == expected_metadata
                @test !haskey(model.metadata, "epoch")
                @test !haskey(model.metadata, "datum")
                metadata_payload = synthetic_payload()
                metadata_payload["epoch"] = "synthetic frozen epoch"
                metadata_payload["datum"] = "synthetic test reference"
                annotated = construct_from_payload(directory, metadata_payload)
                @test annotated.metadata["epoch"] == "synthetic frozen epoch"
                @test annotated.metadata["datum"] == "synthetic test reference"
            end

            @testset "Independent analytic interpolation and units" begin
                for (h, a, b) in (
                    (125, 20.0, 40.0), (675.0, -30.0, 145.0),
                    (500.0, 0.0, 250.0), (0.0, -60.0, 0.0),
                    (1000.0, 70.0, 300.0), (250.0, -10.0, 80.0),
                    (0.0, 5.0, 130.0), (1000.0, 5.0, 130.0),
                )
                    actual = GridAPI.density_state(model, h, deg2rad(a), deg2rad(b))
                    check_state(actual, analytic_fields(h, a, b))
                end
            end

            @testset "Cyclic longitude seam" begin
                h, a = 125.0, 20.0
                left = analytic_fields(h, a, 300.0)
                right = analytic_fields(h, a, 0.0)
                expected = ntuple(i -> (left[i] + right[i]) / 2, 5)
                for longitude in (330.0, -30.0, 1050.0, -750.0)
                    check_state(GridAPI.density_state(model, h, deg2rad(a), deg2rad(longitude)), expected)
                end
                for longitude in (0.0, 360.0, -360.0, 1080.0)
                    check_state(GridAPI.density_state(model, h, deg2rad(a), deg2rad(longitude)),
                                analytic_fields(h, a, 0.0))
                end
            end

            @testset "Frozen time and stored east north up wind" begin
                expected = analytic_fields(125.0, 20.0, 40.0)
                for elapsed in (-12, 0.0, 3600.0), wind in (false, true)
                    check_state(GridAPI.density_state(model, 125, deg2rad(20.0), deg2rad(40.0), elapsed, wind), expected)
                end
                @test any(!iszero, Tuple(GridAPI.density_state(model, 125, 0.0, 0.0, 0, false)[3]))
            end

            @testset "Explicit domain policy and existing tolerance" begin
                for (h, a) in ((-1.0, 0.0), (1001.0, 0.0), (500.0, -60.001), (500.0, 70.001))
                    @test_throws DomainError GridAPI.density_state(model, h, deg2rad(a), 0.0)
                end
                # The existing checked segment uses 1e-12 times the axis scale.
                check_state(GridAPI.density_state(model, 1000.0 + 5e-10, 0.0, 0.0), analytic_fields(1000.0, 0.0, 0.0))
                check_state(GridAPI.density_state(model, -5e-10, 0.0, 0.0), analytic_fields(0.0, 0.0, 0.0))
                check_state(GridAPI.density_state(model, 500.0, deg2rad(70.0) + 5e-13, 0.0), analytic_fields(500.0, 70.0, 0.0))
                @test_throws DomainError GridAPI.density_state(model, 1000.0 + 1e-8, 0.0, 0.0)
                @test_throws DomainError GridAPI.density_state(model, -1e-8, 0.0, 0.0)

                vacuum = GridAPI.GRAMGridAtmosphereModel(;
                    planet="mars", surrogate_file=source, above_grid=:vacuum, vacuum_temperature=175,
                )
                check_state(GridAPI.density_state(vacuum, 1001.0, 0.0, 0.0), (0.0, 175.0, 0.0, 0.0, 0.0))
                check_state(GridAPI.density_state(vacuum, 1000.0, 0.0, 0.0), analytic_fields(1000.0, 0.0, 0.0))
                @test_throws DomainError GridAPI.density_state(vacuum, -1.0, 0.0, 0.0)
                @test_throws DomainError GridAPI.density_state(vacuum, 1001.0, deg2rad(70.1), 0.0)
                @test_throws DomainError GridAPI.density_state(vacuum, 1001.0, deg2rad(-60.1), 0.0)
                @test_throws ArgumentError GridAPI.GRAMGridAtmosphereModel(; planet="mars", surrogate_file=source, above_grid=:clamp)
                for bad_temperature in (0.0, -1.0, NaN, Inf)
                    @test_throws ArgumentError GridAPI.GRAMGridAtmosphereModel(; planet="mars", surrogate_file=source, vacuum_temperature=bad_temperature)
                end
                for component in 1:4, bad in (NaN, Inf, -Inf)
                    coordinates = [125.0, 0.0, 0.0, 0.0]
                    coordinates[component] = bad
                    @test_throws DomainError GridAPI.density_state(model, coordinates...)
                    @test_throws DomainError GridAPI.density_state(vacuum, coordinates...)
                end
            end

            @testset "Deep-copy, serialization and fresh file loading" begin
                expected = analytic_fields(125.0, 20.0, 40.0)
                duplicate = deepcopy(model)
                @test duplicate.surrogate.rho !== model.surrogate.rho
                @test duplicate.source_sha256 == actual_sha
                check_state(GridAPI.density_state(duplicate, 125.0, deg2rad(20.0), deg2rad(40.0)), expected)
                io = IOBuffer()
                serialize(io, model)
                seekstart(io)
                restored = deserialize(io)
                @test restored.source_sha256 == actual_sha
                @test restored.metadata == model.metadata
                check_state(GridAPI.density_state(restored, 125.0, deg2rad(20.0), deg2rad(40.0)), expected)

                changed_payload = synthetic_payload()
                changed_payload["fields"]["density_kgm3"] .= 2.0
                write_payload(source, changed_payload)
                reloaded = GridAPI.GRAMGridAtmosphereModel(; planet="mars", surrogate_file=source)
                @test reloaded.source_sha256 != actual_sha
                @test GridAPI.density_state(reloaded, 125.0, 0.0, 0.0)[1] == 2.0
                check_state(GridAPI.density_state(model, 125.0, deg2rad(20.0), deg2rad(40.0)), expected)
                write_payload(source)
            end
        end

        @testset "Resolver precedence, missing files, LFS and SHA verification" begin
            first_root = mktempdir(directory)
            second_root = mktempdir(directory)
            first_file = joinpath(first_root, "mars_surrogate.jls")
            second_file = write_payload(joinpath(second_root, "mars_surrogate.jls"))
            @test GridAPI.resolve_gram_grid_file("mars"; search_roots=[first_root, second_root]) == second_file
            write_payload(first_file)
            @test GridAPI.resolve_gram_grid_file("mars"; search_roots=[first_root, second_root]) == first_file
            @test GridAPI.resolve_gram_grid_file("mars"; surrogate_file=source, search_roots=[first_root, second_root]) == source
            @test_throws ArgumentError GridAPI.resolve_gram_grid_file("mars"; surrogate_file=joinpath(directory, "missing.jls"), search_roots=[second_root])
            @test_throws ArgumentError GridAPI.resolve_gram_grid_file("pluto"; surrogate_file=source)

            if Sys.iswindows()
                # Windows symlinks may require privileges unavailable in CI.
                @test_skip !Sys.iswindows()
            else
                rm(first_file)
                symlink(joinpath(first_root, "missing-target.jls"), first_file)
                @test_throws ArgumentError GridAPI.resolve_gram_grid_file("mars"; search_roots=[first_root, second_root])
                rm(first_file)
            end

            pointer = "version https://git-lfs.github.com/spec/v1\noid sha256:$(repeat("a", 64))\nsize 12345\n"
            write(first_file, pointer)
            message = argument_error_text(() -> GridAPI.resolve_gram_grid_file("mars"; search_roots=[first_root, second_root]))
            @test occursin(r"(?i)lfs|pointer", message)
            @test occursin(first_file, message)
            message = argument_error_text(() -> GridAPI.GRAMGridAtmosphereModel(; planet="mars", surrogate_file=first_file, search_roots=[second_root]))
            @test occursin(r"(?i)lfs|pointer", message)
            @test_throws ArgumentError GridAPI.GRAMGridAtmosphereModel(; planet="mars", surrogate_file=source, expected_sha256=repeat("0", 64))
            @test_throws ArgumentError GridAPI.GRAMGridAtmosphereModel(; planet="mars", surrogate_file=source, expected_sha256="invalid")
            @test_throws ArgumentError GridAPI.GRAMGridAtmosphereModel(; planet="venus", surrogate_file=source)
        end

        @testset "Reject invalid grid schema and values" begin
            # Valid serialized objects that do not implement the schema.
            for payload in (42, [], Dict{String,Any}())
                argument_error_text(() -> construct_from_payload(directory, payload))
            end
            for key in ("status", "type", "planet", "grid", "fields")
                payload = synthetic_payload()
                delete!(payload, key)
                argument_error_text(() -> construct_from_payload(directory, payload))
            end
            for (key, value) in (("status", "error"), ("type", "other"), ("planet", "venus"), ("grid", nothing), ("fields", nothing))
                payload = synthetic_payload()
                payload[key] = value
                argument_error_text(() -> construct_from_payload(directory, payload))
            end
            for (axis, values) in (
                ("alt_km", Float64[]), ("alt_km", [0.0]),
                ("alt_km", [0.0, 0.0, 1.0]), ("alt_km", [1.0, 0.25, 0.0]),
                ("alt_km", [0.0, NaN, 1.0]), ("alt_km", [0.0, 0.25, Inf]),
                ("lat_deg", [-91.0, -10.0, 70.0]), ("lat_deg", [-60.0, -10.0, 91.0]),
                ("lat_deg", [-60.0, -10.0, -10.0]), ("lat_deg", [-60.0, NaN, 70.0]),
                ("lon_deg", [-1.0, 80.0, 210.0, 300.0]),
                ("lon_deg", [0.0, 80.0, 210.0, 360.0]),
                ("lon_deg", [0.0, 80.0, 80.0, 300.0]),
                ("lon_deg", [0.0, 210.0, 80.0, 300.0]),
                ("lon_deg", [0.0, 80.0, 210.0, Inf]),
            )
                payload = synthetic_payload()
                payload["grid"][axis] = values
                argument_error_text(() -> construct_from_payload(directory, payload))
            end
            for key in ("alt_km", "lat_deg", "lon_deg")
                payload = synthetic_payload()
                delete!(payload["grid"], key)
                argument_error_text(() -> construct_from_payload(directory, payload))
            end
            for key in SYNTHETIC_FIELD_KEYS
                payload = synthetic_payload()
                delete!(payload["fields"], key)
                argument_error_text(() -> construct_from_payload(directory, payload))
                for bad in (NaN, Inf)
                    payload = synthetic_payload()
                    payload["fields"][key][1] = bad
                    argument_error_text(() -> construct_from_payload(directory, payload))
                end
                payload = synthetic_payload()
                payload["fields"][key] = ones(2, 3, 4)
                argument_error_text(() -> construct_from_payload(directory, payload))
            end
            for replacement in (ones(3, 3), ones(3, 3, 4, 1), "invalid")
                payload = synthetic_payload()
                payload["fields"]["density_kgm3"] = replacement
                argument_error_text(() -> construct_from_payload(directory, payload))
            end
            for (key, bad) in (("density_kgm3", -0.1), ("temperature_K", 0.0), ("temperature_K", -1.0))
                payload = synthetic_payload()
                payload["fields"][key][1] = bad
                argument_error_text(() -> construct_from_payload(directory, payload))
            end
            zero_density = synthetic_payload()
            zero_density["fields"]["density_kgm3"] .= 0.0
            @test GridAPI.density_state(construct_from_payload(directory, zero_density), 125.0, 0.0, 0.0)[1] == 0.0
            for contents in (UInt8[], UInt8[0xff, 0x00, 0x23, 0x14], Vector{UInt8}(codeunits("not a grid")))
                path = joinpath(directory, "invalid_binary.jls")
                write(path, contents)
                message = argument_error_text(() -> GridAPI.GRAMGridAtmosphereModel(; planet="mars", surrogate_file=path))
                @test occursin(path, message)
            end
        end

        @test GridAPI._GRAM_WRAPPER[] === native_wrapper_before
        @test GridAPI._GRAM_WRAPPER_FILE[] == native_wrapper_file_before
    end
end

end # module GridAtmosphereTests
