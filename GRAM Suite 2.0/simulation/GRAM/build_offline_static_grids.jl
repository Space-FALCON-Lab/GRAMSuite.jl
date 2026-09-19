#!/usr/bin/env julia
# Compatibility entry point. All generator implementation lives in scripts/.
let legacy_root = normpath(dirname(dirname(@__DIR__)))
    canonical = normpath(joinpath(legacy_root, "..", "scripts", "build_offline_static_grids.jl"))
    if haskey(ENV, "SPACEAGORA_GRAM_ROOT")
        Base.include(@__MODULE__, canonical)
    else
        # Capture the old layout during include without changing the caller's
        # environment. The canonical builder stores this resolved root.
        withenv("SPACEAGORA_GRAM_ROOT" => legacy_root) do
            Base.include(@__MODULE__, canonical)
        end
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
