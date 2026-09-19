# Generate a frozen atmosphere from a recorded recipe

Most SpaceAGORA users should load a surrogate and run their guidance, control or
autonomy algorithm. Generating a new surrogate is an advanced workflow that
requires a separate native GRAM installation. Installing this package and
validating a recipe do not load or install native GRAM.

The generator samples one planet at one frozen instant. It stores density,
temperature and nominal east/north/up winds on a spatial grid. A simulation's
elapsed time does not advance this atmosphere. Use native GRAM when the study
requires atmospheric time variation, or generate and validate another frozen
snapshot for the scenario being studied.

## Start with a native-free plan

From the GRAMSuite.jl checkout, with an existing parent directory for run records:

```sh
julia --startup-file=no scripts/generate_grid.jl \
  --recipe=recipes/odyssey_p20_frozen_v1.toml \
  --output=/path/to/runs/odyssey-plan
```

The default mode is `dry-run`. It checks the complete recipe, retains its exact
bytes and input lock, records generator/wrapper source identities and writes the
normalized builder arguments. It does not need the native library or planetary
data and does not make a grid. The output directory must be new.

The supplied recipe records the previously evaluated Odyssey P20 diagnostic fine
configuration: 100 to 260 km ellipsoidal height, 40 to 90 degrees geodetic north,
full periodic east-positive longitude, and 1 km by 1 degree by 2.5 degree spacing.
The grid has 161 by 51 by 144 nodes and Float64 fields. Its frozen UTC instant is
2001-11-07T11:51:04.794789Z, in planet-event frame 0 with elapsed time zero.
Mars radii, MOLA off, MapYear 2 followed by the min-max 1 parameter push, seed
1001 and disabled perturbations are explicit. Dust settings retain the pinned
native defaults. Unprobed settings remain labelled as such.

This is a diagnostic recipe, not a released scientific preset. The reference
grid checksum identifies the historical artifact; it is not an expected checksum
of a new serialized file, which also contains a generation timestamp. Independent
checks must establish configuration fidelity and compare the generated values.

## Check the supplied native inputs without loading them

```sh
julia --startup-file=no scripts/generate_grid.jl \
  --recipe=recipes/odyssey_p20_frozen_v1.toml \
  --mode=preflight \
  --gram-root="/path/to/GRAM Suite 2.0" \
  --lib="/path/to/libGRAM.dylib" \
  --spice="/path/to/selected/SPICE" \
  --output=/path/to/runs/odyssey-preflight
```

Use the library filename for your platform. The external GRAM root must contain
`Julia/GRAM.jl` and the selected planet data in `Mars/data` or `Earth/data` when
applicable. Earth additionally requires an explicit `--earth-merra2-path` with
`All Mean/MERRA2All_MM.bin` for the recipe's month. No alternate dataset is chosen
when an explicit path is wrong.

Preflight checks files and hashes; it does not call the native binding. The
Odyssey recipe pins its accepted native-library SHA256 and a sidecar with the
nine retained Mars input identities and thirty supplied SPICE kernel identities.
It rejects extra, missing or changed files in those selected inventories. The
selected data and SPICE directories must match those recorded sets, rather than
arbitrary larger installation directories containing additional files. The old
evidence recorded these selected files, not every auxiliary file in the original
installation. An isolated external root can link the selected files into the
expected layout; the run must still pass native configuration and field checks.

The example's library identity is from the accepted macOS run. Another platform,
native build or input set needs a separately identified recipe and its own
configuration checks. Removing a pin creates a run that captures observed
inputs; it does not establish equivalence to the locked Odyssey run.

## Explicit native generation

After arranging exclusive use of the native installation, repeat the preflight
command with `--mode=native` and a new output directory. This is the only mode
that starts a builder subprocess. It fixes the native library path in both the
builder arguments and the child's `GRAM_LIB`, passes `SPACEAGORA_GRAM_ROOT`, and
uses one Julia/OpenMP/BLAS thread. The child works inside the private run directory
so relative native diagnostics are retained there. No native assets are copied
into the package.

The child receives an explicit `--project` pointing at a newly recorded empty
project and `JULIA_LOAD_PATH=@:@stdlib` (the platform equivalent on Windows). The
generator uses Julia standard libraries, and the retained GRAM binding is
self-contained. The parent's active project, `JULIA_PROJECT` and user package
environment are not inherited as dependency sources. Project bytes are hashed;
no manifest or dependency installation is produced. A different binding that
requires external Julia packages needs a separately reviewed environment change.

Each completed record includes the original recipe and lock, normalized
invocation, Julia command prefix and explicit environment overrides, source
hashes and available Git revision, binding hashes, the native-library hash,
selected planet-data/SPICE inventories, builder log, applied binding arguments,
reported native model version when available, output checksums and execution
status. The reported model version is separate from the recipe's declared label.
Source and asset identities are checked again around execution. Inventories
describe supplied directories; they are not a claim that the native model opened
every file or loaded every listed kernel. No raw environment, credentials,
native source or data binaries are included in the provenance record.

The runner checks axes, field types, finite values, positive density/temperature
and recorded binding settings before promoting the run directory. It rejects
overwrites. A failed run is retained as a separate `.failed-<id>` directory and
the requested successful output location remains absent. `SHA256SUMS` covers
every retained file except itself. A completed run has status
`generated_not_scientifically_accepted`; successful execution does not set
scientific acceptance limits or publish a preset.

## Further planets and validation

Schema 1 supports explicit recipes for Earth, Mars, Venus, Titan, Jupiter,
Uranus and Neptune. It currently uses the uniform Float64 generation path,
full periodic longitude, a UTC frozen instant and nominal winds. Mars requires
explicit radii and settings. For other planets the builder's native radii/default
model settings must be declared as defaults with their version and unprobed
status; this does not make them validated planetary presets. No settings or
coverage for the next planet are supplied implicitly.

Choose the next scenario's coverage and forcing deliberately. Compare a new grid
with independently configured native GRAM at held-out locations, including
longitude seams, poles and altitude transitions. Keep interpolation error and
frozen-time error separate. Check their effects on forces, heating and trajectory
outputs before deciding whether the preset meets that application's needs.
Odyssey exact-pole winds retain longitude-labelled local ENU samples, with no
unique physical-vector claim, averaging or newly invented cutoff.

The lower-level `scripts/build_offline_static_grids.jl` retains direct and
resampling workflows, including nonuniform/adaptive options. Recipe schema 1
does not silently invoke those options. Resampling preserves source atmospheric
metadata and rejects unsupported source bounds; it cannot change the source's
epoch or forcing by changing generation flags.

## Software checks

```sh
julia --startup-file=no test/generator_builder_tests.jl
julia --startup-file=no test/generation_recipe_tests.jl
```

These tests use stand-ins, including a real Julia subprocess that emits synthetic
data. They test validation, input pins, provenance, failure retention and the
builder contract without executing native GRAM. They provide software evidence,
not new atmospheric or trajectory validation.
