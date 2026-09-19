#!/usr/bin/env julia

include(joinpath(@__DIR__, "grid_generation_recipe.jl"))
if abspath(PROGRAM_FILE) == @__FILE__
    GridGenerationRecipe.main(ARGS)
end
