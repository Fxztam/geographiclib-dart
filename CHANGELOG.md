# Changelog

## 2.0.0 - Triaxial Ellipsoid

Date: 2026-04-24

### Added

- Added the Dart triaxial geodesic implementation under `lib/src/triaxial/`.
- Added the public triaxial barrel export `lib/triaxial.dart`.
- Added dataset-driven triaxial validation via `tool/geod3test.dart`.
- Added `test/triaxial_test.dart` for unit, regression, inverse, and direct triaxial checks.
- Added the authoritative reference dataset file `test/data/Geod3Test-v1.txt` required for triaxial validation.
- Added transfer and validation documentation in `README_triaxial.md` and `README_triaxial_transfer.md`.

### Validation

- `dart test test\triaxial_test.dart`
- `dart run tool\geod3test.dart --quick`

Validated on the source integration branch as well:

- full `Geod3Test-v1` run over `500000` cases
- inverse `max = 168.00 eps`
- inverse `mean = 1.258 eps`
- `alp1_max = 36.081 peps`
- `alp2_max = 38.099 peps`
- direct `max = 1674.44 eps`

### Notes

- This release introduces the triaxial ellipsoid implementation and its reference validation workflow.
- The target repository transfer was validated with unit/regression tests and the quick dataset harness before GitHub handoff.

## 1.0.0

- Initial release.
- Dart port of Charles Karney's GeographicLib (geodesic, DMS, Transverse Mercator, UTM).
- 271 tests passing.
