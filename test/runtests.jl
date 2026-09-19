using Test, Libdl
import GRAMSuite

const NATIVE_LIBRARIES_BEFORE = Set(Libdl.dllist())
@testset "GRAMSuite native-free package suite" begin
@test GRAMSuite._GRAM_WRAPPER[] === nothing
include("grid_atmosphere_tests.jl")
include("grid_contract_tests.jl")
include("cache_regression.jl")
include("constructor_recipe_tests.jl")
include("mapyear_profile_tests.jl")
@test GRAMSuite._GRAM_WRAPPER[] === nothing
@test isempty(filter(p -> occursin(r"(?i)(libgram|libcspice)",p),
    collect(setdiff(Set(Libdl.dllist()),NATIVE_LIBRARIES_BEFORE))))
println("All tests use synthetic grids or explicit Julia stand-ins; native GRAM was not loaded.")
end
