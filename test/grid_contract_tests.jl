module GridContractTests
using Test, Serialization, StaticArrays
import GRAMSuite
using ..GridAtmosphereTests: synthetic_payload, write_payload, analytic_fields

@testset "Full-grid format and concurrent owned snapshots" begin
    mktempdir() do directory
        payload=synthetic_payload()
        payload["type"]="full_grid"
        payload["format"]="spaceagora_gram_static_grid_v1"
        path=write_payload(joinpath(directory,"full_grid.jls"),payload)
        model=GRAMSuite.GRAMGridAtmosphereModel(planet="mars",surrogate_file=path)
        @test model.metadata["type"]=="full_grid"
        @test model.metadata["format"]=="spaceagora_gram_static_grid_v1"
        expected=analytic_fields(125.0,20.0,40.0)
        actual=GRAMSuite.density_state(model,125.0,deg2rad(20.0),deg2rad(40.0))
        @test all(isapprox.(Float64[actual[1],actual[2],actual[3]...],collect(expected);rtol=2e-13,atol=2e-13))
        # Full-grid support is deliberately limited to the native-free constructor.
        @test_throws ArgumentError GRAMSuite._gram_load_offline_surrogate(path,"mars")
        for format in (nothing,"unknown_format","spaceagora_gram_surrogate_trilinear_v1")
            invalid=deepcopy(payload)
            format===nothing ? delete!(invalid,"format") : (invalid["format"]=format)
            bad=write_payload(joinpath(directory,"bad_grid.jls"),invalid)
            @test_throws ArgumentError GRAMSuite.GRAMGridAtmosphereModel(planet="mars",surrogate_file=bad)
        end
        legacy=synthetic_payload()
        legacy["format"]="spaceagora_gram_surrogate_trilinear_v1"
        legacy_path=write_payload(joinpath(directory,"surrogate.jls"),legacy)
        @test GRAMSuite.GRAMGridAtmosphereModel(planet="mars",surrogate_file=legacy_path).metadata["format"]==legacy["format"]
        legacy["format"]="spaceagora_gram_static_grid_v1"
        write_payload(legacy_path,legacy)
        @test_throws ArgumentError GRAMSuite.GRAMGridAtmosphereModel(planet="mars",surrogate_file=legacy_path)
        arrays_before=deepcopy((model.surrogate.rho,model.surrogate.T,model.surrogate.wind_e,
            model.surrogate.wind_n,model.surrogate.wind_u))
        outputs=Vector{Any}(undef,256)
        Threads.@threads for i in eachindex(outputs)
            outputs[i]=GRAMSuite.density_state(model,125.0,deg2rad(20.0),deg2rad(40.0),Float64(i),isodd(i))
        end
        @test all(==(actual),outputs)
        @test arrays_before==(model.surrogate.rho,model.surrogate.T,model.surrogate.wind_e,
            model.surrogate.wind_n,model.surrogate.wind_u)
        @test GRAMSuite._GRAM_WRAPPER[]===nothing
        # Serialization contains the snapshot; evaluation does not reopen its source.
        buffer=IOBuffer();serialize(buffer,model);seekstart(buffer)
        rm(path)
        restored=deserialize(buffer)
        @test GRAMSuite.density_state(restored,125.0,deg2rad(20.0),deg2rad(40.0))==actual
        duplicate=deepcopy(restored)
        duplicate.surrogate.rho[1]=123.0
        @test restored.surrogate.rho[1]!=123.0
        @test GRAMSuite.density_state(restored,125.0,deg2rad(20.0),deg2rad(40.0))==actual
    end
end
end # module GridContractTests
