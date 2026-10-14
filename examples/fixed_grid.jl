# A small API demonstration. These constants are invented, not NASA GRAM output.
using GRAMSuite
using Serialization
using SHA

mktempdir() do directory
    payload = Dict(
        "status" => "ok",
        "type" => "surrogate_trilinear",
        "planet" => "mars",
        "provenance" => "Synthetic constant fields for the getting-started example",
        "grid" => Dict(
            "alt_km" => [0.0, 1.0],
            "lat_deg" => [-45.0, 45.0],
            "lon_deg" => [0.0, 120.0, 240.0],
        ),
        "fields" => Dict(
            "density_kgm3" => fill(0.02, 2, 2, 3),
            "temperature_K" => fill(210.0, 2, 2, 3),
            "wind_ew_ms" => fill(5.0, 2, 2, 3),
            "wind_ns_ms" => zeros(2, 2, 3),
            "wind_up_ms" => zeros(2, 2, 3),
        ),
    )
    file = joinpath(directory, "synthetic_mars.jls")
    serialize(file, payload)
    digest = bytes2hex(open(sha256, file))
    model = GRAMGridAtmosphereModel(
        planet="mars", surrogate_file=file, expected_sha256=digest,
    )
    rho, temperature, wind_enu = density_state(
        model, 500.0, deg2rad(10.0), deg2rad(60.0), 0.0, true,
    )
    @assert isapprox(rho, 0.02)
    @assert isapprox(temperature, 210.0)
    @assert all(isapprox.(wind_enu, (5.0, 0.0, 0.0)))
    println("Density: ", rho, " kg/m^3")
    println("Temperature: ", temperature, " K")
    println("East/north/up wind: ", Tuple(wind_enu), " m/s")
    println("Synthetic example complete; no native GRAM or atmosphere data was used.")
end
