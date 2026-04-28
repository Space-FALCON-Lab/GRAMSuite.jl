# GRAMSuite.jl

A Julia package for querying planetary atmosphere models using NASA's Global Reference Atmosphere Model (GRAM) Suite. Supports atmosphere density, temperature, and wind queries for Earth, Mars, Venus, Titan, Jupiter, Uranus, and Neptune.

## Prerequisites

- Julia 1.10 or later
- NASA GRAM Suite 2.0 (see [Acquiring GRAM](#acquiring-gram))
- A C build toolchain:
  - **macOS:** GNU Make (`brew install make`)
  - **Linux:** `make` and a C compiler (GCC)
  - **Windows:** MSYS2/MinGW with `mingw32-make`

## Acquiring GRAM

GRAM Suite 2.0 is controlled software distributed by NASA. Request a copy through the NASA Software Catalog:

**https://software.nasa.gov/software/MFS-33888-1**

Once approved, you will receive an archive containing the GRAM Suite source code, data files, and SPICE kernels. Place the extracted `GRAM Suite 2.0` folder at the root of this repository (alongside `src/` and `Project.toml`).

## Setup

### 1. Build the shared library

From the repository root, run the provided build script and pass the path to your `GRAM Suite 2.0` folder:

**macOS / Linux:**
```bash
cd "GRAM Suite 2.0/simulation/GRAM"
./build_gram.sh "$(pwd)/../.."
```

**Windows (PowerShell):**
```powershell
cd "GRAM Suite 2.0\simulation\GRAM"
.\build_gram.ps1 (Resolve-Path "..\..").Path
```

**Windows (CMD):**
```cmd
cd "GRAM Suite 2.0\simulation\GRAM"
build_gram.cmd "%CD%\..\.."
```

The script detects your platform, runs `setup_cspice.sh`, and produces the shared library:

| Platform | Output |
|----------|--------|
| Linux    | `GRAM Suite 2.0/Build/lib/libGRAM.so` |
| macOS    | `GRAM Suite 2.0/Build/lib/libGRAM.dylib` |
| Windows  | `GRAM Suite 2.0/Build/lib/libGRAM.dll` |

It also writes a `gram.env` file (shell) with `GRAM_ROOT` and `GRAM_LIB` set for the current machine.

### 2. Install the Julia package

See [SpaceAGORA.jl documentation](https://www.spacefalconlab.com/SpaceAGORA.jl/) for package-specific installation.

From the repository root in the Julia REPL:

```julia
using Pkg
Pkg.develop(path=".")
```

Or add it to another project:

```julia
Pkg.add(url="https://github.com/your-org/GRAMSuite.jl")
```

### 3. Verify the installation

Run the smoke test (requires `GRAM_ROOT` to be set or the library already built):

```bash
cd "GRAM Suite 2.0/simulation/GRAM"
./run_julia_smoke_test.sh
```

Expected output includes `GRAM smoke test passed` and sample atmospheric state values.

## Usage

```julia
using GRAMSuite

# Construct a model for a planet.
# gram_directory and gram_data_directory default to "GRAM Suite 2.0" relative
# to the repository root; override with absolute paths as needed.
model = GRAMAtmosphereModel(
    planet_name = "mars",
    initial_time = InitialTime(year=2024, month=6, day=1)
)

# Query atmosphere state at a point.
# h: altitude in metres, lat/lon in radians, elapsed_time in seconds
h_m       = 50_000.0
lat_rad   = deg2rad(22.0)
lon_rad   = deg2rad(48.0)
elapsed_s = 0.0

rho, T, wind = density_state(model, h_m, lat_rad, lon_rad, elapsed_s, true)
# rho  — density (kg/m³)
# T    — temperature (K)
# wind — SVector{3} [east, north, up] wind components (m/s)
```

### Path overrides

If your GRAM Suite is not at the default location, pass explicit paths:

```julia
model = GRAMAtmosphereModel(
    planet_name        = "earth",
    gram_root_directory = "/path/to/GRAM Suite 2.0",
    gram_data_directory = "/path/to/GRAM Suite 2.0",
    gram_library_path   = "/path/to/libGRAM.so"   # optional
)
```

Alternatively, set environment variables before starting Julia:

```bash
export GRAM_ROOT="/path/to/GRAM Suite 2.0"
export GRAM_LIB="/path/to/GRAM Suite 2.0/Build/lib/libGRAM.so"
```

## Configuration

Runtime behaviour is controlled through environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `GRAM_ROOT` | auto-detected | Path to the `GRAM Suite 2.0` root directory |
| `GRAM_LIB` | `<GRAM_ROOT>/Build/lib/libGRAM.*` | Explicit path to the shared library |
| `SPACEAGORA_GRAM_WIND_MODE` | `auto` | Wind field to return: `nominal` (mean) or `perturbed` (stochastic) |
| `SPACEAGORA_GRAM_GLOBAL_LOCK` | `on` | Thread-safety lock for GRAM calls; set to `off` only for single-threaded use |
| `SPACEAGORA_GRAM_OFFLINE_SURROGATE` | `off` | Use pre-built surrogate payloads instead of live GRAM calls: `on`, `off`, or `auto` |
| `SPACEAGORA_GRAM_STATIC_GRID` | `false` | Use a frozen in-memory interpolation grid instead of live GRAM calls |

## SPICE kernels

The GRAM Suite ships with a starter pack of NAIF SPICE kernels covering dates from 2000-01-01 to 2100-01-01. If your mission requires a different date range or higher-precision ephemeris data, download the appropriate kernels from the [NAIF website](https://naif.jpl.nasa.gov/naif/data.html) and point GRAM to them via its SPICE configuration.

## License

This wrapper is released under the MIT License. NASA GRAM Suite 2.0 is subject to its own export-controlled distribution terms; refer to the documentation included with your GRAM distribution.
