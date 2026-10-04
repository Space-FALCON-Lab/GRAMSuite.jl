# Compare horizontal wind interpolation with reference samples

Use `scripts/validation/horizontal_wind.py` to compare a rectangular grid with disjoint reference observations at the same ellipsoid height. This diagnostic requires Python 3.9 or newer and its standard library. It does not load GRAM, generate native observations, change the atmosphere backend or select a release tolerance.

A surrogate can have above-ground interpolation nodes and still miss horizontal wind structure. This tool helps distinguish the effect of horizontal spacing from the effect of explicitly chosen vector-frame treatments. It reports every requested treatment without automatically selecting the best one.

## Prepare two CSV files

Generate and freeze the node and reference point plans separately before collecting reference values. Use the same atmospheric configuration for both files, and retain its recipe, native/model/data identities, epoch, coordinate conventions and collector revision alongside the observations. CSV hashes in the output identify the observations but do not establish that their generation settings are correct.

Both CSVs require these columns:

| Column | Meaning |
| --- | --- |
| `id` | Nonempty sample identifier, unique within its file |
| `slice_id` | Identifier shared by one node plane and its reference samples |
| `altitude_km` | Geodetic ellipsoid height in kilometres, identical within a slice |
| `latitude_deg` | Geodetic latitude in degrees |
| `longitude_deg` | East-positive longitude in degrees, in `[0,360)` |
| `wind_east_ms`, `wind_north_ms`, `wind_up_ms` | Stored local wind components in m/s |
| `clearance_km` | Independently established native surface clearance in kilometres |
| `planetocentric_latitude_deg` | Required only for the planetocentric Cartesian treatment; use the actual native readback |

Each node slice must form a complete rectangular latitude/longitude grid. Nonuniform axes are supported. Each reference sample must be inside its matching node bounds, at exactly the same height, and disjoint from all generation-node coordinates. Renaming a node does not make it a held-out sample. The tool rejects duplicate rows/IDs, incomplete grids, nonfinite fields, mismatched slices/heights, overlapping validation coordinates and out-of-domain queries. This first utility supports bounded longitude tiles narrower than 180 degrees and excludes exact poles. It does not validate periodic seams or the native pole-wind convention.

The files may contain multiple slices at different heights or locations. There is no vertical interpolation between them. Stored components must have consistent units and conventions; the tool cannot recover a native wind basis from unlabeled numbers.

## Run a comparison

```sh
python3 scripts/validation/horizontal_wind.py \
  --nodes nodes.csv --reference reference.csv --out comparison \
  --methods stored_enu cartesian_geodetic cartesian_planetocentric \
  --strides 1 2 4 8
```

For a 9-by-9 node slice, strides 1, 2, 4 and 8 use 9-by-9, 5-by-5, 3-by-3 and 2-by-2 grids respectively. Each stride must retain both endpoints of each axis. Omitting `--strides` compares the full grid. Omitting `--methods` evaluates only `stored_enu`.

- `stored_enu` interpolates the three stored components independently.
- `cartesian_geodetic` rotates each node's vector into Cartesian coordinates using a geodetic local basis, interpolates, then rotates into the query's geodetic local basis.
- `cartesian_planetocentric` follows the same procedure using the explicitly supplied planetocentric latitudes.

The two Cartesian choices are diagnostic hypotheses. Neither establishes the native model's wind basis or should silently replace the production convention. A lower error for one choice in a finite sample does not prove that choice is physically correct.

The output directory must be new. Validation failures occur before it is created. `errors.csv` contains signed component errors, vector norms and query/corner clearances. `summary.json` records input/tool SHA256 identities and max, linearly interpolated 95th-percentile and RMS errors by slice and treatment. It reports the full cohort and the subset with strictly positive clearance at the query and all positively weighted corners. No weight cutoff is used; zero-weight corners are ignored. A clearance estimate is not a conservative terrain uncertainty bound.

Compare cohort counts as well as errors. Different strides may include different above-ground subsets, so their subset maxima need not describe the same sample population. Keep the full-cohort comparisons available. Large errors, worsening on refinement, or failed hypotheses are results to retain rather than reasons to discard observations.

The tool always records `accuracy_requirements: null` and `release_accepted: false`. Application requirements belong in a prospectively approved validation protocol. If these references influence a subsequent grid design or method choice, treat them as development data and collect fresh acceptance points for the frozen candidate. Direct agreement with native GRAM is a model-reconstruction result, not independent physical validation.

## Native-free checks

```sh
python3 test/horizontal_wind_tests.py
```

The tests cover analytic ENU and Cartesian fields, nonuniform grids, coordinate transforms, disjoint validation, clearance cohorts, malformed inputs and complete CLI output. They run in CI without native assets or third-party Python packages. Julia runtime and generator tests continue separately.
