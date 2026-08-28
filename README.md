# GRAMSuite.jl

A Julia package for querying planetary atmosphere models using NASA's Global Reference Atmosphere Model (GRAM) Suite. Supports atmosphere density, temperature, and wind queries for Earth, Mars, Venus, Titan, Jupiter, Uranus, and Neptune.

## Prerequisites

- Julia 1.10 or later
- NASA GRAM Suite 2.0 (see [Acquiring GRAM](#acquiring-gram))
- Git LFS (`git lfs install`) — this repository stores its `.jls` surrogate
  payloads and SPICE kernels via LFS
- A C/C++ toolchain and GNU Make:
  - **Linux:** `gcc`/`g++` and `make`
  - **macOS:** Clang from the Xcode command-line tools, plus GNU Make
    (`brew install make` — the build files use GNU Make features that the
    `/usr/bin/make` shipped with macOS does not support, so `gmake` is required)
  - **Windows:** MSYS2/MinGW with `mingw32-make`
- CSPICE, which is normally already handled for you:
  - **Linux x86_64 and Windows MinGW** — an archive ships inside the GRAM Suite
    distribution and is selected automatically. Nothing to install.
  - **macOS** — `brew install cspice`
  - **Linux arm64/aarch64** — no bundled archive; see
    [CSPICE on unbundled platforms](#cspice-on-unbundled-platforms)

You do **not** need a Fortran compiler. GRAM's build configuration names
`gfortran`, but `Build/makefile.defs` ends with `undefine FC`, which disables
the Fortran example targets — the shared-library build never invokes it.

## Acquiring GRAM

GRAM Suite 2.0 is controlled software distributed by NASA. Request a copy through the NASA Software Catalog:

**https://software.nasa.gov/software/MFS-33888-1**

Once approved, you will receive an archive containing the GRAM Suite source code, data files, and SPICE kernels. Place the extracted `GRAM Suite 2.0` folder at the root of this repository (alongside `src/` and `Project.toml`).

## Setup

### 1. Build the shared library

The Julia package calls GRAM through FFI, so the native shared library must be
compiled once per machine. It is never shipped prebuilt, and a library built on
one machine is not usable on another.

One command does the whole thing. The path argument is optional — the script
walks up from its own location to find the tree containing `Build/` and
`Julia/`:

**macOS / Linux:**
```bash
"GRAM Suite 2.0/simulation/GRAM/build_gram.sh"
```

**Windows (PowerShell):**
```powershell
& ".\GRAM Suite 2.0\simulation\GRAM\build_gram.ps1"
```

**Windows (CMD):**
```cmd
"GRAM Suite 2.0\simulation\GRAM\build_gram.cmd"
```

Pass an explicit root only if this folder has been moved outside a GRAM Suite
tree — quote it, since the default folder name contains a space:
`./build_gram.sh "/absolute/path/to/GRAM Suite 2.0"`.

The script, in order:

1. runs `Build/setup_cspice.sh` to put the correct CSPICE archive in place —
   **do not run this yourself**
2. runs `make shared` with one job per core, using the right make binary for
   the host (`make`, `gmake`, or `mingw32-make`)
3. produces the shared library:

   | Platform | Output |
   |----------|--------|
   | Linux    | `GRAM Suite 2.0/Build/lib/libGRAM.so` |
   | macOS    | `GRAM Suite 2.0/Build/lib/libGRAM.dylib` |
   | Windows  | `GRAM Suite 2.0/Build/lib/libGRAM.dll` |

4. writes `simulation/GRAM/gram.env` (`gram.env.ps1` on Windows) exporting
   `GRAM_ROOT` and `GRAM_LIB` for this machine
5. writes `simulation/GRAM/.gram-build-manifest`, recording the host platform
   and root path the artifacts were built for

A build from clean takes roughly 20 seconds on 24 cores, or about 3 minutes of
total CPU time.

`gram.env` and `.gram-build-manifest` contain absolute paths for the machine
that produced them. They are gitignored and must never be committed — a
committed one makes the next person's build fail against a path that does not
exist on their machine.

#### Building by hand

The wrapper is the supported path. If you need to drive `make` yourself:

```bash
cd "GRAM Suite 2.0/Build"
./setup_cspice.sh
make clean && make shared -j        # gmake on macOS, mingw32-make on Windows
```

`make shared` is the only target that produces `libGRAM`. Plain `make` builds
the static per-planet libraries and the example executables but **not** the
shared library, and neither do the per-planet targets — `make Mars -j` leaves
`Build/lib` with `libMars.a` and no `libGRAM.so`. There is no single-planet
shortcut to a smaller `libGRAM`; the `shared` target links every planet.

Driving `make` by hand also writes no `.gram-build-manifest`, so the staleness
detection described below has nothing to compare against on the next wrapper
run.

#### Rebuilding

Add `--clean` (`-Clean` in PowerShell) to force a full rebuild:

```bash
"GRAM Suite 2.0/simulation/GRAM/build_gram.sh" --clean
```

You usually do not need it. The build manifest records the host tag and root
path, and the script cleans automatically when either has changed — so moving
or renaming the tree on the same machine is handled for you, and so is a tree
that arrived from a different operating system.

Reach for `--clean` when a tree was copied from another machine **with
`Build/lib` already populated**, since an existing library of the right name
can otherwise be mistaken for a native one.

#### CSPICE on unbundled platforms

`setup_cspice.sh` normalizes the archive name the makefiles expect:

- **Linux x86_64** — bundled `common/cspice/lib/cspice_gcc85.a` is copied to
  `cspice_linux_x86_64.a`
- **Windows MinGW** — bundled `common/cspice/lib/cspice_mingw64.a` is used as-is
- **macOS** — taken from Homebrew (`brew install cspice`)
- **Linux arm64/aarch64** — no bundled archive. Install CSPICE so that one of
  `/usr/lib/libcspice.a`, `/usr/local/lib/libcspice.a`, `/usr/lib64/libcspice.a`,
  or `/usr/lib/aarch64-linux-gnu/libcspice.a` exists, or build with an explicit
  override:

  ```bash
  cd "GRAM Suite 2.0/Build"
  make shared -j SPICE_LIB=/absolute/path/to/cspice.a
  ```

#### Troubleshooting the build

| Symptom | Cause and fix |
|---|---|
| `Build finished but shared library not found`, naming a path that is not yours | A `gram.env`/`.gram-build-manifest` from another machine is present. Delete both and re-run. |
| `Could not auto-find a Linux CSPICE archive` | Architecture with no bundled archive — see above. |
| `gmake not found` on macOS | `brew install make`. GNU Make is required. |
| `make` succeeded but there is no `libGRAM` | You ran `make` or a per-planet target instead of `make shared`. |
| Library builds, but a planet errors at query time | A data problem, not a build problem. Run `git lfs pull` and confirm the planet's `data/` folder came across from the GRAM distribution. |

### 2. Install the Julia package

See [SpaceAGORA.jl documentation](https://www.spacefalconlab.com/SpaceAGORA.jl/) for package-specific installation.

From the repository root in the Julia REPL:

```julia
using Pkg
Pkg.develop(path=".")
```

Or add it to another project:

```julia
Pkg.add(url="https://github.com/Space-FALCON-Lab/GRAMSuite.jl")
```

### 3. Verify the installation

Run the smoke test. It sources the `gram.env` written by the build step, so
run the build first:

```bash
"GRAM Suite 2.0/simulation/GRAM/run_julia_smoke_test.sh"
```

To run Julia directly instead, source that file yourself:

```bash
source "GRAM Suite 2.0/simulation/GRAM/gram.env"
julia "GRAM Suite 2.0/simulation/GRAM/julia_smoke_test.jl"
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

## Height and latitude conventions

The wrapper's `density_state` path accepts SI inputs — height in meters, latitude and
longitude in radians — and converts them once at the native boundary (km / degrees).
The conventions those inputs are interpreted under matter at the kilometer level:

- **Latitude/height pair**: callers supply Bowring **planetodetic** latitude with
  **geodetic height above the reference ellipsoid** (what SpaceAGORA's `rtolatlong`
  produces). The wrapper passes `is_planetocentric=false`, so native GRAM converts
  the pair to planetocentric internally (`Position::convertToPlanetocentric`).
  Labeling geodetic inputs planetocentric mis-references the height by the local
  ellipsoid geometry — up to ~2 km at Mars polar latitudes, which is a ~25% density
  error at aerobraking heights (density scale height ~7 km).
- **Mars height reference**: native MarsGRAM defaults to `isMolaHeights=true`
  (heights above the MOLA areoid, planetocentric inputs required). The wrapper
  defaults this to **false** so heights are referenced to the ellipsoid, matching
  the geodetic inputs above. Opting back in via `mars_mola_heights=true` requires
  planetocentric positions and is therefore incompatible with the wrapper's
  `density_state` inputs — native GRAM raises an error at the first query.
- **Sanity check**: native GRAM echoes the resolved radius via `get_position`
  (`totalRadius`); with the defaults above a known planetocentric radius round-trips
  to sub-meter accuracy. Passing radius-based altitudes (`r − R_equatorial`) as the
  height is wrong under *either* setting — at 80° latitude it lands ~17.5 km below
  the intended point (~14× density).

## License

This wrapper is released under the MIT License. NASA GRAM Suite 2.0 is subject to its own export-controlled distribution terms; refer to the documentation included with your GRAM distribution.
