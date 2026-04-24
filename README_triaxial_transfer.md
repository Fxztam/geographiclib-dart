# Triaxial Transfer Checklist

Date: 2026-04-24

## Purpose

This is the short transfer checklist for the Dart triaxial geodesic files.

Use it immediately before copying or importing the triaxial implementation into another repository.

## Working Directory

Run all commands from:

```powershell
cd D:\works\work-flutter\geographiclib
```

## Reference Dataset

The repository should not carry the large authoritative dataset file.

- Origin: https://doi.org/10.5281/zenodo.12510796
- Dataset version: `1.0`
- Upstream file name: `Geod3Test-v1.txt`
- Expected local placement for the harness: `test/data/Geod3Test-v1.txt`

If the dataset is stored somewhere else locally, run the harness with `--data PATH` or set `GEOD3TEST_DATA`.

## Required Dart Commands

### 1. Unit and Regression Test

```powershell
dart test test\triaxial_test.dart
```

Expected:

- all tests pass

### 2. Quick Dataset Smoke Test

```powershell
dart run tool\geod3test.dart --quick
```

Expected:

- inverse and direct both pass

### 3. Full Dataset Gate

```powershell
dart run tool\geod3test.dart --report-every 50000
```

Expected:

- full `Geod3Test-v1` run passes
- inverse max error stays `<= 500 eps`
- inverse mean error stays `<= 20 eps`
- direct max error stays `<= 5000 eps`

## Strongly Recommended Extra Checks

### Focused inverse anchor

```powershell
dart run tool\geod3test.dart --start 5728 --count 1 --inv
```

### Focused direct anchor

```powershell
dart run tool\geod3test.dart --start 7928 --count 3 --dir
```

### Larger slice

```powershell
dart run tool\geod3test.dart --start 0 --count 50000 --report-every 5000
```

## Current Verified Full-Run Result

Latest validated full-run result:

- `500000` cases
- inverse max: `168.00 eps`
- inverse mean: `1.258 eps`
- `alp1_max = 36.081 peps`
- `alp2_max = 38.099 peps`
- direct max: `1674.44 eps`

## Transfer Rule

Do not transfer the Dart triaxial files unless at least these three commands are green:

1. `dart test test\triaxial_test.dart`
2. `dart run tool\geod3test.dart --quick`
3. `dart run tool\geod3test.dart --report-every 50000`

## Recommendation

Transfer `tool/geod3test.dart` together with the triaxial implementation if possible. It is the most important executable reference check for future validation in the target repository.
