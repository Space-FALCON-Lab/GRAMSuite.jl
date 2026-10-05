# Native-free near-surface Mars atmosphere

`GRAMNearSurfaceAtmosphereModel` evaluates a frozen Mars-GRAM lower atmosphere from 5 m above the local terrain to
the payload's areoid-height top: 75 km in version 1.0.0 of the published preset, and 81 km in version 1.1.0. It
reproduces native Mars-GRAM's own near-surface rule from stored component fields, and returns
density, temperature and pressure. It needs no native GRAM installation, SPICE kernel or GRAM data file.
- Version 1 payloads (format `spaceagora_mars_near_surface_scalars_v1`, preset versions 1.0.0 and 1.1.0) store no winds.
- A version 2 payload (format `spaceagora_mars_near_surface_v2`, planned as preset version 1.2.0) has the same scalars
  plus a wind layer: east, north and vertical winds, and the speed of sound that limits them ([Winds](#winds-format-2)).

```julia
using GRAMSuite

model = GRAMNearSurfaceAtmosphereModel(; surrogate_file="mars_near_surface.jls",
    expected_sha256="<the payload's published SHA256>")

# geodetic latitude and east longitude in degrees, height above the reference ellipsoid in metres
state = near_surface_state(model, -4.5, 137.4, 250.0)
state.density_kgm3, state.temperature_K, state.pressure_Pa, state.regime, state.surface_layer_model_status

# the fixed-grid calling convention: height in metres, angles in radians; the wind vector is the stored winds of a
# version 2 payload, and zero for a version 1 payload
rho, T, wind = density_state(model, 250.0, deg2rad(-4.5), deg2rad(137.4))

# version 2 payloads: the scalars plus the winds and the speed of sound
near_surface_winds_available(model)
w = near_surface_wind_state(model, -4.5, 137.4, 250.0)
w.wind_east_ms, w.wind_north_ms, w.wind_up_ms, w.sound_speed_ms, w.wind_clipped, w.wind_regime
```

The payload must be supplied explicitly; nothing is searched for or downloaded. The published preset
`mars_global_near_surface_p20_frozen_v1` is distributed as GRAMSuite release assets in two versions, 1.0.0 (to 75 km)
and 1.1.0 (to 81 km). Each archive's README states its validation and terms. SpaceAGORA's `surrogate_preset_model`
installs and verifies a named version.

## Evaluation rule

The query is converted to planetocentric latitude and radius on the reference ellipsoid. The terrain component gives the
query's own surface height `zs` and areoid radius, 2D linear on Mars-GRAM's MOLA lattice, as native evaluates them. The
first exposed table level is `L1 = max(floor(zs + 0.3) + 1, -5)` km. It selects one of three regimes:

| Regime | Heights | Rule |
|---|---|---|
| D3 | at or above `L1` | temperature and gas constant linear, and ln p linear, between the bracketing table levels |
| D4 | from 30 m above the surface to `L1` | temperature linear from the 30 m value to the `L1` level value; `p = p_L1 exp((L1 - z) / H)`, with `H = (T_L1 + T_5m) / Q_L` |
| D5 | 5 to 30 m above the surface | temperature linear from the 5 m to the 30 m value; the D4 pressure law |

**Levels.** The table levels are 1 km apart up to 10 km and 5 km apart up to 75 km. Version 1.1.0 adds the upper
table's first two levels, at 80.0323 and 85.0323 km areoid height, where native Mars-GRAM places them at the frozen
instant. The D3 rule applies between them as between any two levels.

**Components.** The level states, the 5 m and 30 m surface temperatures and the surface-layer coefficient `Q_L` are
stored on a lattice and interpolated bilinearly. `Q_L` has a fitted form in each native 7.5-degree latitude band and
9-degree surface cell.

**Edges.** A query on a band or cell edge belongs to both sides; if its own model is missing, the neighbour across that
edge is used. At a corner, where a band edge meets a cell edge, the diagonal neighbour is tried after the two edge
neighbours.

**Status.** Every surface-layer result reports its model's status: `"qualified"`, or a provisional status that the
payload records. At or above `L1` the status is `"not used"`.

## Winds (format 2)

A version 2 payload stores pre-clip east and north winds and sound-speed offsets as evaluated outputs on native-aligned
rows and longitudes. No GRAM table is stored. The scalars are the version 1 rule's, unchanged.

**Height.** Five regimes interpolate linearly in height between stored endpoints:

| Regime | Heights | Endpoints |
|---|---|---|
| D5 | 5 to 30 m above the surface | the 5 m and 30 m surface fields |
| D4 | 30 m above the surface to `L1` | the 30 m field and level `L1` |
| D3 | `L1` to 75 km | the bracketing lower-table levels |
| D2 | 75 km to 80 + ho km | the 75 km level and the upper table's first level, 80 + ho km (80.0323 at the frozen instant) |
| D1 | 80 + ho to 81 km | the upper table's first two levels |

**Horizontal.** Each endpoint field is linear in latitude and longitude between stored nodes.
- The east wind is stored on both sides of 78.75 N and 82.5 N, where native's pole factor jumps. A query exactly on
  either latitude takes the northern side.
- The surface east wind is stored on both sides of the 0 E seam. A query exactly at 0 E takes the western side.

**Slope and vertical winds.** Native's slope-wind term (`MarsWindEndpoints`' `slope_wind`, unchanged) adds the slope
winds and the vertical wind. It uses the terrain slopes from ±0.25 degree differences of the payload's terrain.

**Speed of sound.** c = sqrt(γ P / ρ), with P and ρ from the scalars.
- γ comes from a reference on native's 50 K and pressure-decade cells, plus the stored offset.
- At or below 80.0 km areoid height, water vapour at native's 20 % relative humidity is mixed in.
- Above 80.0 km, a composition term at the query's own molecular weight is added, without water.
- **The 80.0 km switch:** the side is chosen by the query's own areoid height. Native's height arithmetic differs in its
  last bits, so within 1e-9 km of 80.0 km the side may differ from native's. `near_surface_wind_state` flags such
  queries (`in_switch_band`).

**Clipping.** The east and north winds are clipped to ±0.7 c, as native clips them. The vertical wind is not.
`near_surface_wind_state` also returns the unclipped components and whether a component was clipped.

**Refusals.** A query needing a stored node that is unavailable, or a speed of sound that the reference or the
composition does not determine, throws `DomainError` naming the field, endpoint, knot and longitude, or the cell.
Winds are never returned as zero for lack of data.

## Refusals

A query outside the supported domain throws `DomainError` with the first reason that applies:
1. planetocentric latitude beyond the payload's limit (85 degrees for the published preset);
2. surface height at or above its volcano limit (9 km);
3. clearance below its minimum (5 m);
4. areoid height above its top (75 km in 1.0.0, 81 km in 1.1.0);
5. a needed component is unavailable;
6. no surface-layer model exists for the query's cell.

There is no extrapolation and no native fallback.

## Conventions

- The latitude is geodetic, on the payload's reference ellipsoid. The height is measured above that ellipsoid.
- The longitude is east-positive and periodic.
- Query elapsed time must be finite but is not applied: the atmosphere is a frozen snapshot.
- `density_state` returns the stored winds of a version 2 payload, whatever its `wind` argument selects; as with the
  stored winds of the grid presets, callers mask them. For a version 1 payload it returns a zero wind vector, because
  that payload stores no winds.
- The constructor hashes the payload before and after deserialization. It rejects a SHA256 mismatch, a file that changed while loading, a Git LFS pointer, and a payload whose format, planet, array shapes or support limits are invalid.
- Treat the model's arrays and metadata as read-only. Concurrent queries may share one model.

## Data terms

Published near-surface payloads are licensed under CC BY 4.0; see [DATA_LICENSE.md](../DATA_LICENSE.md), including its
note on the terrain component. The code is MIT licensed.
