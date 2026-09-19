# GRAM Simulation Drop-In

This folder is a portable bridge for using GRAM from a simulation project.

## What this gives you

- `build_gram.sh`: Bash build script for macOS/Linux/MSYS shells.
- `run_julia_smoke_test.sh`: Bash Julia smoke test.
- `build_gram.ps1`: native Windows PowerShell build script.
- `run_julia_smoke_test.ps1`: native Windows PowerShell smoke test.
- `build_gram.cmd` and `run_julia_smoke_test.cmd`: Windows CMD wrappers.
- `julia_smoke_test.jl`: the actual Julia smoke test.
- `gram.env` / `gram.env.ps1`: generated env files with `GRAM_ROOT` and `GRAM_LIB`.

## Quick start

macOS/Linux/MSYS shell:

1. `./build_gram.sh /absolute/path/to/GRAM Suite 2.0`
2. `./run_julia_smoke_test.sh`

Windows PowerShell:

1. `.\build_gram.ps1 "C:\path\to\GRAM Suite 2.0"`
2. `.\run_julia_smoke_test.ps1`

Windows CMD:

1. `build_gram.cmd "C:\path\to\GRAM Suite 2.0"`
2. `run_julia_smoke_test.cmd`

## Notes

- macOS uses `gmake` (GNU make). Linux uses `make`.
- Windows uses `mingw32-make` (or `make`) from MSYS2/MinGW.
- Scripts call `Build/setup_cspice.sh` automatically.
- On Windows, bundled `common/cspice/lib/cspice_mingw64.a` is used.

## Offline atmosphere generation

Use the [recipe generation guide](../../../docs/grid_generation.md) for a recorded frozen snapshot. It defines the planet, epoch, coordinates, forcing settings and spatial coverage, and records source/input checksums. Ordinary users load an already prepared surrogate through the [native-free grid API](../../../README.md#native-free-fixed-grids).

The generator has one implementation at `scripts/build_offline_static_grids.jl`. The entry point in this folder forwards to it for existing commands. Keep the repository layout intact, or invoke the canonical script with `SPACEAGORA_GRAM_ROOT` pointing to the external native installation.

From the repository root, inspect the supported advanced options without loading native GRAM:

```sh
julia --startup-file=no scripts/build_offline_static_grids.jl --help
```

Direct generation and resampling remain explicit advanced workflows. Resampling preserves the source atmosphere and rejects unsupported coverage. Adaptive source resampling requires matching altitude endpoints; use fixed bands for a supported subset. No automatic per-planet altitude policy or automatic Earth January override is implemented. Set the epoch, bounds and Earth MERRA2 inputs explicitly in the generation recipe.

The supplied Odyssey recipe is a diagnostic configuration. New planetary presets need their own selected coverage, recorded native configuration and held-out comparisons before being advertised as supported.

## Runtime selection

`GRAMGridAtmosphereModel` is the native-free fixed-grid interface. It uses explicit coverage handling and has no native fallback. `GRAMAtmosphereModel` retains advanced native/hybrid capabilities. The [package README](../../../README.md) describes those interfaces and their distinct configuration requirements; generator success does not choose a simulation backend.

