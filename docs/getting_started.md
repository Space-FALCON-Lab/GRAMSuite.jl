# Getting started with GRAMSuite

GRAMSuite provides atmosphere models that you can query from Julia or use in a
SpaceAGORA simulation. A frozen model reads a saved snapshot and interpolates its
fields. The native backend calls a separately installed NASA GRAM Suite. Start
with the frozen model if you want to learn the interface or use an existing,
validated snapshot.

## Install the Julia package

Follow the [standalone installation](../README.md#installation) with Julia 1.10
or later. If your goal is a SpaceAGORA mission, use Julia 1.12 and the
[SpaceAGORA route below](#use-it-with-spaceagora), which selects the package
revision recorded by that checkout.

Package installation does not generate an atmosphere. A real study also needs
a trusted payload whose epoch, coordinates, coverage and validation are suitable
for the mission. Native execution additionally needs the official GRAM runtime
and its inputs; see the [native setup](../README.md#advanced-native-backend).

## Make a first query

Use the standalone installation above for this example, then run the following
in that Julia environment. The script accompanies this documentation update;
older package revisions, including SpaceAGORA's existing pinned wrapper, may
not contain it. If you are keeping a SpaceAGORA-pinned environment, continue
with the [SpaceAGORA walkthrough](#use-it-with-spaceagora) instead.

```julia
using GRAMSuite
include(joinpath(pkgdir(GRAMSuite), "examples", "fixed_grid.jl"))
```

The [example](../examples/fixed_grid.jl) creates a tiny temporary grid of invented
constant values, loads it, and queries a point at 500 metres altitude. It prints
density near 0.02 kg/m³, temperature near 210 K and east/north/up winds near
(5, 0, 0) m/s. Small floating-point rounding differences are normal. The temporary
file is removed when the example finishes. No atmosphere download, native library
or SPICE kernel is involved. This demonstrates the API, not a physical Mars model.

For an actual snapshot, replace the temporary input with your trusted `.jls`
file and its published checksum:

```julia
model = GRAMGridAtmosphereModel(
    planet="mars",
    surrogate_file="/absolute/path/to/mars_surrogate.jls",
    expected_sha256="<the payload's 64-character SHA256>",
)
rho, temperature, wind_enu = density_state(
    model, 153_000.0, deg2rad(77.5), deg2rad(315.0), 0.0, true,
)
```

The coordinates must lie inside that payload's domain. Height is in metres,
angles in radians and elapsed time in seconds; the result uses kg/m³, K and m/s.
Use the payload's stated height and latitude conventions. Elapsed time does not
advance a frozen atmosphere, and the wind flag does not replace its stored winds.

## Use it with SpaceAGORA

SpaceAGORA keeps GRAMSuite as a separate repository at `data/GRAMSuite.jl` and
connects it through `SpaceAGORAGRAMSuiteExt`. Loading both Julia packages activates
that extension. The folder alone does not install a package, a native library or
an atmosphere payload.

For a first complete mission, follow SpaceAGORA's
[Odyssey surrogate walkthrough](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/src/tutorials/odyssey_surrogate.md).
From the root of a SpaceAGORA Git checkout, its setup and run commands are:

```sh
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate_env/setup.jl
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl
```

The dedicated environment pins GRAMSuite and uses that SpaceAGORA checkout.
The first run retrieves the named frozen atmosphere and the scenario's public
kernels and gravity coefficients, about 224 MB in total; subsequent runs reuse
the installed assets. Read the walkthrough before running it, including its
output-directory, offline and SPICE-session requirements. This is a frozen-grid
workflow; it does not build or start native GRAM.

In that environment, `using SpaceAGORA` followed by `import GRAMSuite` makes
SpaceAGORA's atmosphere adapters available. The
[atmosphere guide](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/src/user/atmosphere_models.md)
shows how to select a named preset or supply a grid to a mission configuration.
Use `import GRAMSuite` here to keep the standalone and SpaceAGORA adapter names
distinct. For native missions, use the separate
[GRAMSuite setup guide](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/src/user/gramsuite_setup.md).

## Choose a model and check compatibility

| Route | Requirements and behavior |
| --- | --- |
| Standalone Julia package | Julia 1.10 or later; repository CI exercises 1.10 and 1.12. `Project.toml` currently identifies package version 0.1.0. Record the resolved Git revision as well. |
| SpaceAGORA integration | Julia 1.12 and the GRAMSuite revision pinned by the selected SpaceAGORA setup. Keep that environment's Project and Manifest files; a package version alone does not identify its source. |
| Fixed grid | `GRAMGridAtmosphereModel` reads supported trilinear or full-grid payloads. Queries stay within the recorded altitude/latitude domain; longitude is periodic. |
| Near-surface Mars | `GRAMNearSurfaceAtmosphereModel` uses a separately supplied near-surface payload. Payload format and preset version determine coverage and whether winds are stored; see the [model guide](near_surface_atmosphere.md). |
| Native backend | `GRAMAtmosphereModel` needs a compatible compiled NASA GRAM Suite installation and planetary inputs. The [native guide](../README.md#advanced-native-backend) covers platform setup. |

The Julia package version, the payload format and the named preset version are
different identifiers. For example, near-surface preset 1.2.0 is not GRAMSuite
package version 1.2.0. Upgrading the Julia code does not change a saved atmosphere's
epoch or extend its validity.

## Validation and attribution

The synthetic example checks installation and the query interface. Repository
tests exercise interpolation, boundaries and data contracts; they do not certify
an atmosphere for a mission. A checksum verifies file identity. It does not prove
physical accuracy. Consult each payload's validation record, and assess both
interpolation error and the effect of freezing the atmosphere in time.

The Julia wrapper is developed by the SpaceAGORA contributors and is
[MIT licensed](../LICENSE). NASA GRAM Suite is the source model for generated
atmosphere products and has its own distribution terms. Published generated
payloads carry the terms in [DATA_LICENSE.md](../DATA_LICENSE.md); cite the preset
name and version, credit NASA GRAM, and retain the terrain or scenario-asset
attribution that applies. The near-surface terrain component also credits MOLA.

For further work, see [grid generation](grid_generation.md),
[terrain geometry](terrain_geometry.md), [wind endpoints](mars_wind_endpoints.md)
and [horizontal-wind diagnostics](horizontal_wind_validation.md).
