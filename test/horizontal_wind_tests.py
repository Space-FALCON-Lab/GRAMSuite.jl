import contextlib
import csv
import importlib.util
import io
import json
import math
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('horizontal_wind', ROOT/'scripts/validation/horizontal_wind.py')
wind = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wind)


def point(identifier, latitude, longitude, cartesian=False):
    vector = wind.transform((10., -20., 30.), latitude, longitude, inverse=True) if cartesian else (latitude+longitude, 2*latitude-longitude, 3.)
    return dict(id=identifier, slice_id='test', altitude_km=5., latitude_deg=latitude, longitude_deg=longitude,
                planetocentric_latitude_deg=latitude, wind_east_ms=vector[0], wind_north_ms=vector[1], wind_up_ms=vector[2], clearance_km=1.)


def fixture(cartesian=False):
    nodes = [point(f'n{i}{j}', i/2, 10+j/2, cartesian) for i in range(3) for j in range(3)]
    reference = [point('a', .23, 10.31, cartesian), point('b', .71, 10.82, cartesian)]
    return nodes, reference


def write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class HorizontalWindTests(unittest.TestCase):
    def test_affine_component_field_on_nested_grids(self):
        n, r = fixture()
        values = wind.compare(n, r, ['stored_enu'], [1, 2])
        self.assertEqual(len(values), 4)
        self.assertLess(max(v['vector_error_ms'] for v in values), 1e-12)

    def test_constant_cartesian_field(self):
        n, r = fixture(True)
        values = wind.compare(n, r, ['cartesian_geodetic', 'cartesian_planetocentric'], [1, 2])
        self.assertLess(max(v['vector_error_ms'] for v in values), 1e-12)

    def test_nonuniform_axes(self):
        n = [point(f'n{i}{j}', p, l) for i, p in enumerate((0., .2, 1.)) for j, l in enumerate((10., 10.7, 11.))]
        values = wind.compare(n, [point('r', .6, 10.8)], ['stored_enu'], [1, 2])
        self.assertLess(max(v['vector_error_ms'] for v in values), 1e-12)

    def test_frame_roundtrip(self):
        for p, l in [(-90, 0), (-35, 41), (0, 360), (90, 190)]:
            vector = (7., -20., 31.)
            out = wind.transform(wind.transform(vector, p, l), p, l, inverse=True)
            self.assertLess(math.hypot(*(a-b for a, b in zip(out, vector))), 1e-12)

    def test_positive_weight_corners_and_clearance_cohorts(self):
        n, r = fixture()
        n[0]['clearance_km'] = -1
        values = wind.compare(n, r, ['stored_enu'], [2])
        self.assertTrue(all(not v['all_above_surface'] for v in values))
        summaries = wind.summarize(values)
        clean = next(s for s in summaries if s['cohort'] != 'all')
        self.assertEqual(clean['count'], 0)
        self.assertIsNone(clean['maximum_ms'])
        value, minimum = wind.interpolate(wind.make_grid(n), n[-1], 'stored_enu')
        self.assertEqual(minimum, 1.)

    def test_reference_must_be_disjoint(self):
        n, _ = fixture()
        with self.assertRaisesRegex(ValueError, 'overlap'):
            wind.compare(n, [dict(n[0], id='different-label')], ['stored_enu'], [1])

    def test_incomplete_grid(self):
        n, r = fixture()
        with self.assertRaisesRegex(ValueError, 'complete rectangular'):
            wind.compare(n[:-1], r, ['stored_enu'], [1])

    def test_invalid_stride_and_method(self):
        n, r = fixture()
        for methods, strides in [(['other'], [1]), (['stored_enu'], [0]), (['stored_enu'], [3]), ([], [1]), (['stored_enu'], [1, 1])]:
            with self.assertRaises(ValueError):
                wind.compare(n, r, methods, strides)

    def test_reference_height_and_bounds(self):
        n, r = fixture()
        for change in ({'altitude_km': 6.}, {'latitude_deg': 1.001}, {'longitude_deg': 9.999}):
            with self.assertRaises(ValueError):
                wind.compare(n, [dict(r[0], **change)], ['stored_enu'], [1])

    def test_slice_sets_and_height_consistency(self):
        n, r = fixture()
        with self.assertRaises(ValueError):
            wind.compare(n, [dict(r[0], slice_id='missing')], ['stored_enu'], [1])
        n[0]['altitude_km'] = 6
        with self.assertRaises(ValueError):
            wind.compare(n, r, ['stored_enu'], [1])

    def test_centric_basis_requires_coordinates(self):
        n, r = fixture()
        del r[0]['planetocentric_latitude_deg']
        with self.assertRaisesRegex(ValueError, 'explicit native latitude'):
            wind.compare(n, r, ['cartesian_planetocentric'], [1])

    def test_csv_rejects_bad_rows(self):
        n, _ = fixture()
        variants = [n+[dict(n[0], id='different')], n+[dict(n[0])], [dict(n[0], wind_east_ms=float('nan'))], [dict(n[0], latitude_deg=90)], [dict(n[0], longitude_deg=360)], [dict(n[0], id='')]]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'nodes.csv'
            for rows in variants:
                write(path, rows)
                with self.assertRaises(ValueError):
                    wind.load_samples(path)
            path.write_text('id,id\nx,x\n')
            with self.assertRaises(ValueError):
                wind.load_samples(path)

    def test_cli_outputs_identities_without_claiming_acceptance(self):
        n, r = fixture()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root/'nodes.csv', n)
            write(root/'reference.csv', r)
            args = ['--nodes', str(root/'nodes.csv'), '--reference', str(root/'reference.csv'), '--out', str(root/'out'), '--strides', '1', '2']
            with contextlib.redirect_stdout(io.StringIO()):
                wind.main(args)
            report = json.loads((root/'out/summary.json').read_text())
            self.assertFalse(report['native_executed'])
            self.assertFalse(report['release_accepted'])
            self.assertIsNone(report['accuracy_requirements'])
            self.assertEqual(len(report['nodes_sha256']), 64)
            self.assertEqual(len(report['summaries']), 4)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                wind.main(args)

    def test_validation_failure_creates_no_outputs(self):
        n, _ = fixture()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root/'nodes.csv', n)
            write(root/'reference.csv', [dict(n[0], id='overlap')])
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                wind.main(['--nodes', str(root/'nodes.csv'), '--reference', str(root/'reference.csv'), '--out', str(root/'out')])
            self.assertFalse((root/'out').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
