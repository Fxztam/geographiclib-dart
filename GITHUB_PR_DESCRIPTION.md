# GitHub PR Description

## Suggested Title

Release 2.0.0: add triaxial ellipsoid support

## PR Body

### Summary

This PR introduces the Dart triaxial ellipsoid implementation into `geographiclib_dart` and prepares the package for the `2.0.0` release line.

The transferred implementation includes the core triaxial solver, public exports, dedicated tests, a dataset-driven validation runner, the reference dataset, and repository-local validation documentation.

### What Changed

- added the triaxial implementation under `lib/src/triaxial/`
- added the public barrel `lib/triaxial.dart`
- added `test/triaxial_test.dart`
- added the dataset runner `tool/geod3test.dart`
- added `test/data/Geod3Test-v1.txt`
- added `README_triaxial.md` and `README_triaxial_transfer.md`

### Validation

Validated in the target repository:

- `dart test test\triaxial_test.dart`
- `dart run tool\geod3test.dart --quick`

Validated on the source integration side before transfer:

- `dart run tool\geod3test.dart --report-every 50000`

Full reference result over `Geod3Test-v1`:

- `500000` test cases
- inverse max error: `168.00 eps`
- inverse mean error: `1.258 eps`
- direct max error: `1674.44 eps`

All values are within the current harness thresholds.

### Release Positioning

This change is intended to be represented in the changelog as:

- `2.0.0 - Triaxial Ellipsoid`

### Outcome

The target repository now contains the full Dart triaxial file set, its validation entry points, and the documentation required to rerun the same checks before GitHub publication.
