#!/usr/bin/env python3
"""Compare horizontal wind interpolants with disjoint reference CSV samples.

Python standard library only. Reads observations; never loads native GRAM.
Frame treatments are explicit hypotheses, not declarations of native wind basis.
"""
import argparse
import bisect
import csv
import hashlib
import json
import math
from pathlib import Path

COORDS = ('altitude_km', 'latitude_deg', 'longitude_deg')
WINDS = ('wind_east_ms', 'wind_north_ms', 'wind_up_ms')
NUMERIC = COORDS + WINDS + ('clearance_km',)
METHODS = ('stored_enu', 'cartesian_geodetic', 'cartesian_planetocentric')


def load_samples(path):
    with Path(path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not set(('id', 'slice_id') + NUMERIC) <= set(fields):
            raise ValueError('CSV needs unique headers and the documented sample columns')
        samples = []
        ids, positions = set(), set()
        for raw in reader:
            if None in raw or any(value is None for value in raw.values()):
                raise ValueError('Malformed CSV row')
            row = dict(raw)
            if not row['id'].strip() or row['id'] in ids or not row['slice_id'].strip():
                raise ValueError('Sample IDs must be nonempty and unique; slice_id is required')
            ids.add(row['id'])
            for key in NUMERIC + (('planetocentric_latitude_deg',) if 'planetocentric_latitude_deg' in fields else ()):
                row[key] = float(row[key])
                if not math.isfinite(row[key]):
                    raise ValueError('Sample fields must be finite')
            if not -90 < row['latitude_deg'] < 90 or not 0 <= row['longitude_deg'] < 360:
                raise ValueError('This bounded diagnostic excludes poles and requires east longitude in [0,360)')
            if 'planetocentric_latitude_deg' in row and not -90 < row['planetocentric_latitude_deg'] < 90:
                raise ValueError('Planetocentric latitude must be strictly between the poles')
            position = (row['slice_id'],) + tuple(row[key] for key in COORDS)
            if position in positions:
                raise ValueError('Duplicate position within a slice')
            positions.add(position)
            samples.append(row)
    if not samples:
        raise ValueError('Empty sample file')
    return samples


def basis(latitude, longitude):
    p, l = math.radians(latitude), math.radians(longitude)
    return ((-math.sin(l), math.cos(l), 0.0),
            (-math.sin(p)*math.cos(l), -math.sin(p)*math.sin(l), math.cos(p)),
            (math.cos(p)*math.cos(l), math.cos(p)*math.sin(l), math.sin(p)))


def transform(vector, latitude, longitude, inverse=False):
    axes = basis(latitude, longitude)
    if inverse:
        return tuple(math.fsum(a*b for a, b in zip(vector, axis)) for axis in axes)
    return tuple(math.fsum(vector[k]*axes[k][j] for k in range(3)) for j in range(3))


def make_grid(rows):
    if len({r['altitude_km'] for r in rows}) != 1:
        raise ValueError('Each slice must have one common ellipsoid height')
    lat = sorted({r['latitude_deg'] for r in rows})
    lon = sorted({r['longitude_deg'] for r in rows})
    if len(lat) < 2 or len(lon) < 2 or len(rows) != len(lat)*len(lon):
        raise ValueError('Nodes must form a complete rectangular grid')
    if lon[-1] - lon[0] >= 180:
        raise ValueError('This diagnostic requires a bounded longitude tile narrower than 180 degrees')
    return lat, lon, {(r['latitude_deg'], r['longitude_deg']): r for r in rows}


def bracket(axis, value):
    if not axis[0] <= value <= axis[-1]:
        raise ValueError('Reference sample outside node bounds; extrapolation is unsupported')
    index = min(bisect.bisect_right(axis, value)-1, len(axis)-2)
    return index, (value-axis[index])/(axis[index+1]-axis[index])


def interpolate(grid, query, method):
    lat, lon, nodes = grid
    i, a = bracket(lat, query['latitude_deg'])
    j, b = bracket(lon, query['longitude_deg'])
    weighted = [(nodes[lat[i], lon[j]], (1-a)*(1-b)),
                (nodes[lat[i+1], lon[j]], a*(1-b)),
                (nodes[lat[i], lon[j+1]], (1-a)*b),
                (nodes[lat[i+1], lon[j+1]], a*b)]
    weighted = [(row, weight) for row, weight in weighted if weight > 0]
    lat_key = 'planetocentric_latitude_deg' if method == 'cartesian_planetocentric' else 'latitude_deg'
    transformed = []
    for row, weight in weighted:
        vector = tuple(row[k] for k in WINDS)
        if method != 'stored_enu':
            vector = transform(vector, row[lat_key], row['longitude_deg'])
        transformed.append((vector, weight))
    value = tuple(math.fsum(w*v[k] for v, w in transformed) for k in range(3))
    if method != 'stored_enu':
        value = transform(value, query[lat_key], query['longitude_deg'], inverse=True)
    return value, min(row['clearance_km'] for row, _ in weighted)


def compare(nodes, references, methods, strides):
    if not methods or len(set(methods)) != len(methods) or any(m not in METHODS for m in methods):
        raise ValueError('Choose distinct documented frame treatments')
    if not strides or len(set(strides)) != len(strides) or any(type(s) is not int or s < 1 for s in strides):
        raise ValueError('Node strides must be distinct positive integers')
    if 'cartesian_planetocentric' in methods and any('planetocentric_latitude_deg' not in r for r in nodes+references):
        raise ValueError('Planetocentric treatment needs explicit native latitude at every node and reference')
    source_slices = {r['slice_id'] for r in nodes}
    if source_slices != {r['slice_id'] for r in references}:
        raise ValueError('Node and reference slice sets must match')
    # Overlap is rejected even when sample labels differ.
    positions = {tuple(r[k] for k in COORDS) for r in nodes}
    if any(tuple(r[k] for k in COORDS) in positions for r in references):
        raise ValueError('Reference coordinates overlap generation nodes; use disjoint validation data')
    results = []
    for name in sorted(source_slices):
        full = [r for r in nodes if r['slice_id'] == name]
        truth = [r for r in references if r['slice_id'] == name]
        lat, lon, mapping = make_grid(full)
        if any(r['altitude_km'] != full[0]['altitude_km'] for r in truth):
            raise ValueError('Reference height must equal its generation slice height')
        for stride in strides:
            if (len(lat)-1) % stride or (len(lon)-1) % stride:
                raise ValueError('Stride must retain both endpoints of both axes')
            grid = lat[::stride], lon[::stride], mapping
            for method in methods:
                for row in truth:
                    value, clearance = interpolate(grid, row, method)
                    difference = tuple(value[k]-row[key] for k, key in enumerate(WINDS))
                    magnitude = math.hypot(*difference)
                    if not all(math.isfinite(x) for x in difference) or not math.isfinite(magnitude):
                        raise ValueError('Wind comparison overflowed; check sample magnitudes')
                    results.append(dict(id=row['id'], slice_id=name, node_stride=stride, method=method,
                                        east_error_ms=difference[0], north_error_ms=difference[1], up_error_ms=difference[2],
                                        vector_error_ms=magnitude,
                                        query_clearance_km=row['clearance_km'], minimum_corner_clearance_km=clearance,
                                        all_above_surface=row['clearance_km'] > 0 and clearance > 0))
    return results


def summarize(results):
    groups = sorted({(r['slice_id'], r['node_stride'], r['method']) for r in results})
    summaries = []
    for name, stride, method in groups:
        selected = [r for r in results if (r['slice_id'], r['node_stride'], r['method']) == (name, stride, method)]
        for cohort in ('all', 'above_surface_query_and_corners'):
            values = sorted(r['vector_error_ms'] for r in selected if cohort == 'all' or r['all_above_surface'])
            index = .95*(len(values)-1)
            i = max(0, int(index))
            p95 = (values[i] + (index-i)*(values[min(i+1, len(values)-1)]-values[i])) if values else None
            summaries.append(dict(slice_id=name, node_stride=stride, method=method, cohort=cohort, count=len(values),
                                  maximum_ms=max(values) if values else None, p95_ms=p95,
                                  rms_ms=(max(values)*math.sqrt(math.fsum((v/max(values))**2 for v in values)/len(values)) if max(values) else 0.0) if values else None))
    return summaries


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nodes', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--methods', nargs='+', choices=METHODS, default=['stored_enu'])
    parser.add_argument('--strides', nargs='+', type=int, default=[1])
    args = parser.parse_args(argv)
    if args.out.exists():
        parser.error('Output directory already exists; choose a fresh path')
    try:
        node_bytes, reference_bytes = args.nodes.read_bytes(), args.reference.read_bytes()
        results = compare(load_samples(args.nodes), load_samples(args.reference), args.methods, args.strides)
        if args.nodes.read_bytes() != node_bytes or args.reference.read_bytes() != reference_bytes:
            raise ValueError('Input changed while the comparison was running')
    except (ValueError, OSError, KeyError, ArithmeticError) as error:
        parser.error(str(error))
    args.out.mkdir(parents=True)
    with (args.out/'errors.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    report = dict(schema_version=1, native_executed=False, release_accepted=False, accuracy_requirements=None,
                  native_wind_basis_established=False, methods=args.methods, node_strides=args.strides,
                  nodes_sha256=hashlib.sha256(node_bytes).hexdigest(), reference_sha256=hashlib.sha256(reference_bytes).hexdigest(),
                  tool_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), summaries=summarize(results))
    (args.out/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(f"Compared {len(results)} wind predictions; diagnostic results saved to {args.out}")


if __name__ == '__main__':
    main()
