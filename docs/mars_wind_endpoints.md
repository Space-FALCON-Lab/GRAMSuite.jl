# Bounded frozen Mars wind endpoints

`GRAMSuite.MarsWindEndpoints` evaluates supplied wind observations without native
GRAM or SPICE. It is a separate component: it does not replace grid interpolation,
install a preset, generate data, or change the atmosphere selected by a simulation.

The current supported envelope is planetocentric latitude 33.9 to 34.6 degrees
and east longitude 40 to 42.5 degrees, at one explicitly recorded frozen epoch.
The supplied terrain tile additionally imposes its geodetic/ellipsoid-height
domain. At least 5 m radial terrain clearance and all required physical endpoint
fields are necessary. Unsupported positions and missing fields raise errors.
This local numerical contract does not establish global or mission accuracy.

## Constructing a snapshot

Use an accepted padded `TerrainGeometry.GeometryTile` with declared usable padding
of at least 0.25 degrees in each axis. An `Observation` stores its `FieldID`,
planetocentric latitude, east longitude, areoid height, terrain height and east/north
winds in m/s. All heights are metres. The constructor checks physical exposure and
endpoint placement to 1e-5 m. Snapshot construction also checks observation terrain
against the supplied tile to 1e-4 m. These checks do not recover missing provenance;
the data provider must establish compatible configuration, epoch and coordinates.

`FieldID(:areoid, level_m)` identifies an integer-kilometre areoid field.
`FieldID(:clearance, 5, piece)` and `FieldID(:clearance, 30, piece)` identify
surface-following endpoints. Surface piece 1 is west of and includes 40.5 degrees;
piece 2 is east of it. Observations on the shared boundary must be provided to
both fields with their corresponding identities. No ordinal layer substitution
or below-ground fill is performed.

Each `EndpointField(id, latitudes, charts, observations)` has two fitting latitude
rows. Supply observations in row order, then longitude-node order. Areoid charts
use five harmonic functions; surface charts use ten. `HarmonicChart(lo, hi, n)`
constructs Chebyshev nodes. To preserve recorded generator coordinates, pass
`nodes=recorded_nodes`; these must match the Chebyshev layout within 1e-12 degrees.
Every observation must match its chart node and field identity exactly.

The chart stores a scaled series basis and its actual infinity-norm condition
number, capped at 1000. This condition is not comparable with an unscaled raw
trigonometric matrix. Each chart enforces normalized coordinate magnitude at most
2 and absolute-weight sum at most 10; the combined latitude/longitude weights
also have a cap of 10. Signed weights require every corresponding observation.
Finite testing does not certify continuous coverage for arbitrary supplied data.

Construct the snapshot from the prepared tile and endpoint fields:

```julia
using GRAMSuite.MarsWindEndpoints

snapshot = WindSnapshot(tile, fields;
    solar_offset_hours=recorded_solar_offset,
    epoch=recorded_frozen_epoch,
    time_convention=recorded_time_convention,
    recipe_sha256=verified_recipe_sha256)
```

Epoch and time convention are recorded provenance strings, not parsed calendar
inputs. The recipe checksum identifies the provider's generation recipe; this
constructor does not download or verify that recipe. Changing time requires a
new compatible snapshot. Both surface endpoints and both surface pieces are
required. Missing areoid levels raise an error when a query needs them.

## Querying and limiting winds

```julia
result = winds(snapshot, geodetic_latitude_deg, east_longitude_deg,
               ellipsoid_height_m; sound_speed_m_s=sound_speed)
background = winds(snapshot, geodetic_latitude_deg, east_longitude_deg,
                   ellipsoid_height_m; sound_speed_m_s=sound_speed,
                   include_slope=false)
```

Sound speed must be a separately supplied positive finite value in m/s. Supply
`sound_speed_m_s=nothing` explicitly for unclamped diagnostic output. There is no
native lookup or guessed sound speed. The caller must ensure the supplied sound
speed uses the same atmospheric configuration and query as the wind snapshot.

The query's own terrain determines the first lower level and vertical weights.
Background winds are reconstructed first. Slope winds and vertical wind are then
computed using the unclamped horizontal wind. Finally each horizontal component
is limited to ±0.7 times sound speed. The vertical component is preserved through
that limit. The result exposes `wind_m_s`, `unclamped_m_s`, `background_m_s`,
clamp status, maximum required-field amplification, condition, physical fields,
geometry and provenance.

All stored fields, nodes, coefficients, observations and geometry are immutable.
Caller arrays are copied into immutable tuples at construction. Concurrent reads
share no mutable cache. This component makes no allocation or speedup promise.

Numerical acceptance uses a 1e-8 m/s gate. Float64 radius subtraction and placement
near the 5/30 m boundaries can contribute about 1e-10 m/s; machine-specific final
bits are not a physical effect. Do not infer mission accuracy from this numerical
gate. Other latitude brackets, dates, seams and poles require additional evidence.
