# Native-free near-surface Mars atmosphere

`GRAMNearSurfaceAtmosphereModel` evaluates a frozen Mars-GRAM lower atmosphere from 5 m above the local terrain to
75 km areoid height. It reproduces native Mars-GRAM's own near-surface rule from stored component fields, and returns
density, temperature and pressure. It needs no native GRAM installation, SPICE kernel or GRAM data file. Winds are not
provided.

```julia
using GRAMSuite

model = GRAMNearSurfaceAtmosphereModel(; surrogate_file="mars_near_surface.jls",
    expected_sha256="<the payload's published SHA256>")

# geodetic latitude and east longitude in degrees, height above the reference ellipsoid in metres
state = near_surface_state(model, -4.5, 137.4, 250.0)
state.density_kgm3, state.temperature_K, state.pressure_Pa, state.regime, state.surface_layer_model_status

# the fixed-grid calling convention: height in metres, angles in radians; the wind vector is always zero
rho, T, wind = density_state(model, 250.0, deg2rad(-4.5), deg2rad(137.4))
```

The payload must be supplied explicitly; nothing is searched for or downloaded. The published preset
`mars_global_near_surface_p20_frozen_v1` is distributed as a GRAMSuite release asset. SpaceAGORA's
`surrogate_preset_model` installs and verifies it by name.

## Evaluation rule

The query is converted to planetocentric latitude and radius on the reference ellipsoid. The terrain component gives the
query's own surface height `zs` and areoid radius, 2D linear on Mars-GRAM's MOLA lattice, as native evaluates them. The
first exposed table level is `L1 = max(floor(zs + 0.3) + 1, -5)` km. It selects one of three regimes:

| Regime | Heights | Rule |
|---|---|---|
| D3 | at or above `L1` | temperature and gas constant linear, and ln p linear, between the bracketing table levels |
| D4 | from 30 m above the surface to `L1` | temperature linear from the 30 m value to the `L1` level value; `p = p_L1 exp((L1 - z) / H)`, with `H = (T_L1 + T_5m) / Q_L` |
| D5 | 5 to 30 m above the surface | temperature linear from the 5 m to the 30 m value; the D4 pressure law |

**Components.** The level states, the 5 m and 30 m surface temperatures and the surface-layer coefficient `Q_L` are
stored on a lattice and interpolated bilinearly. `Q_L` has a fitted form in each native 7.5-degree latitude band and
9-degree surface cell.

**Edges.** A query on a band or cell edge belongs to both sides; if its own model is missing, the neighbour across that
edge is used.

**Status.** Every surface-layer result reports its model's status: `"qualified"`, or a provisional status that the
payload records. At or above `L1` the status is `"not used"`.

## Refusals

A query outside the supported domain throws `DomainError` with the first reason that applies:
1. planetocentric latitude beyond the payload's limit (85 degrees for the published preset);
2. surface height at or above its volcano limit (9 km);
3. clearance below its minimum (5 m);
4. areoid height above its top (75 km);
5. a needed component is unavailable;
6. no surface-layer model exists for the query's cell.

There is no extrapolation and no native fallback.

## Conventions

- The latitude is geodetic, on the payload's reference ellipsoid. The height is measured above that ellipsoid.
- The longitude is east-positive and periodic.
- Query elapsed time must be finite but is not applied: the atmosphere is a frozen snapshot.
- `density_state` returns a zero wind vector because the payload stores no winds. Callers that need winds must take them from another source.
- The constructor hashes the payload before and after deserialization. It rejects a SHA256 mismatch, a file that changed while loading, a Git LFS pointer, and a payload whose format, planet, array shapes or support limits are invalid.
- Treat the model's arrays and metadata as read-only. Concurrent queries may share one model.

## Data terms

Published near-surface payloads are licensed under CC BY 4.0; see [DATA_LICENSE.md](../DATA_LICENSE.md), including its
note on the terrain component. The code is MIT licensed.
