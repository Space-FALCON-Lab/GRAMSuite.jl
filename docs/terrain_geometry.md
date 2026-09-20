# Bounded terrain geometry

`GRAMSuite.TerrainGeometry` provides geometry queries without native GRAM,
SPICE, or atmospheric data. Supply a terrain tile and an explicit served domain.
It computes areoid height and radial terrain clearance from a geodetic input.
It does not decide whether an atmosphere model supports that position.

```julia
using GRAMSuite.TerrainGeometry

ellipsoid = Ellipsoid(3_396_190, 3_376_200) # metres
domain = ServedDomain((34, 35), (40, 42.5), (-5_000, 41_250))
latitude = [33.25, 33.75, 34.25, 34.75, 35.25] # planetocentric degrees
longitude = collect(39.488:0.5:42.988)         # east degrees

# Synthetic demonstration fields, not measured Mars terrain.
areoid = fill(3_390_000.0, length(latitude), length(longitude))
terrain = zeros(length(latitude), length(longitude))
tile = GeometryTile(latitude, longitude, areoid, terrain, ellipsoid, domain;
                    minimum_padding_deg=(0.25, 0.25))
result = geometry(tile, 34.5, 41.0, 0.0)
certificate = padding_certificate(tile)
```

The two matrices use latitude as the first dimension and longitude as the second.
Areoid radius and terrain height above the areoid are in metres. Query latitude is
geodetic in degrees; query height is metres above the supplied oblate ellipsoid.
The result reports radius, planetocentric latitude, interpolated areoid radius,
terrain height, areoid height and radial clearance. Radial clearance is radius
minus areoid radius minus terrain height; it is not geodetic terrain altitude.

The constructor snapshots all four arrays into immutable tuple-backed arrays.
Changing caller inputs does not affect a tile. Tile arrays, their views and their
backing storage cannot be changed with array mutation operations. To prepare
different data, use `copy(tile.terrain_height_m)` and construct a new validated
tile. Copies are ordinary editable arrays. A tile can be shared by concurrent
readers without locks. Tuple storage is intended for bounded tiles; this change
does not establish memory or performance targets for global grids.

Served latitude bounds are within [-89, 89] degrees and served ellipsoid height
must be at least -100,000 m. The coordinate helper additionally enforces its
mathematical normal-coordinate chart for the supplied ellipsoid. Tile axes must
increase strictly, exclude poles, and lie in a single longitude chart within
[0, 360). Seam crossing and polar caps are unsupported. Queries outside the
served domain fail; there is no implicit wrapping, clamping or extrapolation.

Padding is checked in planetocentric latitude and east longitude, both in
degrees. The declared positive padding must meet a 64-ulp numerical floor at
each axis magnitude. The constructor subtracts that floor from each measured
margin and rounds down one representable step before comparing to the declared
padding. The certificate exposes raw and guarded margins in south, north, west,
east order, and requirements and floors in latitude, longitude order.

The guard is a bounded numerical policy, not a universal guarantee for every
ellipsoid or math library. It is not terrain uncertainty, a required physical
clearance, or an atmosphere-validity margin. No terrain dataset is bundled with
this component. Existing atmospheric queries and preset retrieval are unchanged.
