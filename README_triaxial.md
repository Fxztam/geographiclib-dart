# Triaxial Test README

Date: 2026-04-24

## Purpose

This document describes the tests that are required for the Dart triaxial geodesic implementation before the triaxial files are transferred into another repository.

The goal is not only to check local unit tests, but to verify the Dart implementation against the authoritative `Geod3Test-v1` reference dataset with the same workflow that was used to validate the JavaScript implementation.

## Triaxial Features

The triaxial module is exposed separately from the standard `Geodesic` API and follows Karney's triaxial naming.

| Class / API | Description |
| --- | --- |
| `Ellipsoid3` | Triaxial ellipsoid parameterization, axis-derived constants, coordinate normalization, and biaxial/prolate/oblate detection |
| `Geodesic3` | Triaxial geodesic solver for the inverse and direct problems |
| `GeodesicLine3` | Triaxial geodesic line representation for stepping along a solved path |
| `Angle` | Internal angle helper used by the triaxial solver |
| `EllipticFunction3` | Elliptic integral and function support required by the triaxial algorithms |
| `triaxial.dart` | Public barrel export for the triaxial API |

The triaxial implementation was added in version `2.0.0`.

## Relevant Dart Files

Core triaxial implementation:

- `lib/src/triaxial/geodesic3.dart`
- `lib/src/triaxial/geodesic_line3.dart`
- `lib/src/triaxial/trigfun.dart`
- `lib/src/triaxial/ellipsoid3.dart`
- `lib/src/triaxial/angle.dart`
- `lib/src/triaxial/elliptic_function3.dart`

Test and validation entry points:

- `test/triaxial_test.dart`
- `tool/geod3test.dart`
- `test/data/Geod3Test-v1.txt`

## Test Types

### 1. Unit and Regression Tests

File:

- `test/triaxial_test.dart`

What this covers:

- basic `Ellipsoid3` construction checks
- basic `Angle` behavior
- `EllipticFunction3` numerical sanity checks
- inverse reference checks against a curated set of dataset cases
- direct reference checks using Cartesian endpoint and direction comparison
- regression anchors for the previously fixed triaxial solver issues

This is the first required test because it is the cheapest executable validation and catches regressions in the ported solver fixes quickly.

Run from the `dart` directory:

```powershell
dart test test\triaxial_test.dart
```

Expected status:

- all tests pass

## 2. Dataset Harness Tests

File:

- `tool/geod3test.dart`

What this does:

- reads `Geod3Test-v1.txt`
- runs inverse validation against reference distance and azimuth values
- runs direct validation against reference endpoint position
- uses Cartesian endpoint comparison for the direct path, which is the correct near-pole metric

Supported flags:

- `--quick`
- `--start N`
- `--count N`
- `--report-every N`
- `--inv`
- `--dir`

Run from the `dart` directory:

```powershell
dart run tool\geod3test.dart --quick
```

## Required Validation Sequence

The following sequence is recommended before transferring the triaxial files into another repository.

### Step 1. Unit and Regression Suite

```powershell
dart test test\triaxial_test.dart
```

Why:

- fastest broad regression check on the ported triaxial solver

### Step 2. Focused Inverse Anchor

```powershell
dart run tool\geod3test.dart --start 5728 --count 1 --inv
```

Why:

- checks the former `A.b` equatorial inverse failure after the `ArcPos0` default branch fix

### Step 3. Focused Direct Anchor

```powershell
dart run tool\geod3test.dart --start 7928 --count 3 --dir
```

Why:

- checks the previously suspicious near-pole direct family with the corrected Cartesian direct metric

### Step 4. Quick Mixed Slice

```powershell
dart run tool\geod3test.dart --quick
```

Why:

- gives a cheap combined inverse/direct smoke test over the first 1000 dataset rows

### Step 5. Larger Reference Slice

```powershell
dart run tool\geod3test.dart --start 0 --count 50000 --report-every 5000
```

Why:

- verifies that the implementation stays stable over a materially larger part of the authoritative dataset

### Step 6. Full Reference Run

```powershell
dart run tool\geod3test.dart --report-every 50000
```

Why:

- this is the final release or transfer gate for the Dart triaxial implementation
- it verifies the full inverse/direct run on all `500000` dataset cases

## Pass Criteria

Current harness thresholds:

- inverse max error: `<= 500 eps`
- inverse mean error: `<= 20 eps`
- direct max error: `<= 5000 eps`

Current fully validated Dart full-run result:

- `500000` cases
- inverse max: `168.00 eps`
- inverse mean: `1.258 eps`
- `alp1_max = 36.081 peps`
- `alp2_max = 38.099 peps`
- direct max: `1674.44 eps`

These values are within the current pass thresholds.

## Working Directory

All Dart commands above are intended to be run from:

```powershell
cd D:\works\work-flutter\geographiclib
```

## Minimal Transfer Gate

Before copying the triaxial Dart files into another repository, the minimum acceptable gate is:

1. `dart test test\triaxial_test.dart`
2. `dart run tool\geod3test.dart --quick`
3. `dart run tool\geod3test.dart --report-every 50000`

If time allows, the focused anchor runs and the 50000-case slice should also be executed because they make failure localization much easier.

## Recommendation

For any future port or repository transfer, keep `tool/geod3test.dart` together with the triaxial files or recreate it immediately in the target repository. Without the dataset harness, it is much harder to distinguish a real solver regression from a local test gap.
