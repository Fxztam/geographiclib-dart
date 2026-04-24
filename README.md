# geographiclib_dart

Dart port of [Charles Karney's GeographicLib](https://geographiclib.sourceforge.io/).

A pure-Dart library — works in any Dart project and in Flutter apps.

## About the port

This is a library in Dart to solve geodesic problems on an ellipsoid model of the Earth with sub-nanometer precision. 

The algorithms are documented in:

> C. F. F. Karney, [Algorithms for geodesics](https://doi.org/10.1007/s00190-012-0578-z), *J. Geodesy* **87**(1), 43–55 (2013); [Addenda](https://geographiclib.sourceforge.io/geod-addenda.html).

It solves the direct and inverse geodesic problems on the WGS84 ellipsoid to better than 15 nm, and provides two independent Transverse Mercator implementations: 6th-order Krüger series ≤ 5 nm and Lee's exact method with full UTM zone handling.

This Dart port covers the complete public API of `geographiclib-geodesic` JS v2.2.0 — supplemented by the Transverse Mercator, UTM, DMS, and EllipticFunction modules from the C++ library. 

Correctness is verified by 271 test cases from Karney's upstream test suites.

Since version `2.0.0`, the package also includes the triaxial ellipsoid implementation as a separate module, following Karney's original structure: the standard ellipsoid-of-revolution API remains under `Geodesic`, while the triaxial solver is exposed separately via `triaxial`.

For the triaxial API, validation workflow, and transfer notes, see [README_triaxial.md](README_triaxial.md).

## Features

| Class / API | Description | Accuracy |
|---|---|---|
| `Geodesic`, `GeodesicLine`, `PolygonArea` | Geodesic direct, inverse, area on WGS84 ellipsoid | < 15 nm |
| `DMS` | DMS string parsing & formatting (Karney API + legacy helpers) | — |
| `KruegerTM`, `TransverseMercatorResult` | Transverse Mercator — 6th-order Krüger series | ≤ 5 nm |
| `TransverseMercatorExact` | Transverse Mercator — Lee's exact method | exact |
| `EllipticFunction` | Jacobi elliptic functions and integrals | — |
| `KarneyUTMConverter`, `KarneyUTMPoint` | UTM ↔ geographic — Karney (≤ 5 nm) | ≤ 5 nm |

**Additional — Snyder formulas** (independent implementation, Karney-free):

| Class / API | Description | Accuracy |
|---|---|---|
| `UTMConverter`, `UTMPoint`, `GeoPoint` | UTM ↔ geographic — Snyder (1987), pp. 57–64 | ~1 mm |

## API parity: geographiclib-geodesic JS v2.2.0 → Dart

All public methods and output flags of the official JavaScript library have been ported 1 : 1 to Dart.

### Geodesic

| JS | Dart |
|---|---|
| `Inverse` | `inverse` |
| `Direct` | `direct` |
| `ArcDirect` | `arcDirect` |
| `Line` | `line` |
| `DirectLine` | `directLine` |
| `ArcDirectLine` | `arcDirectLine` |
| `InverseLine` | `inverseLine` |
| `Polygon` | `polygon` |

### GeodesicLine

| JS | Dart |
|---|---|
| `Position` | `position` |
| `ArcPosition` | `arcPosition` |
| `GenPosition` | `genPosition` |
| `SetDistance` | `setDistance` |
| `SetArc` | `setArc` |
| `GenSetDistance` | `genSetDistance` |

### PolygonArea

| JS | Dart |
|---|---|
| `AddPoint` | `addPoint` |
| `AddEdge` | `addEdge` |
| `Compute` | `compute` |
| `TestPoint` | `testPoint` |
| `TestEdge` | `testEdge` |
| `Clear` | `clear` |

### Output mask flags

| JS | Dart |
|---|---|
| `NONE`, `ALL`, `STANDARD` | `gNone`, `gAll`, `gStandard` |
| `LATITUDE`, `LONGITUDE`, `AZIMUTH` | `gLatitude`, `gLongitude`, `gAzimuth` |
| `DISTANCE`, `DISTANCE_IN` | `gDistance`, `gDistanceIn` |
| `REDUCEDLENGTH`, `GEODESICSCALE` | `gReducedLength`, `gGeodesicScale` |
| `AREA`, `LONG_UNROLL`, `OUT_MASK` | `gArea`, `gLongUnroll`, `gOutMask` |

### Additional functionality from C++ GeographicLib (beyond the JS library)

The Dart port additionally covers algorithms from the C++ GeographicLib that are not part of the JS library:

| Module | Description |
|---|---|
| `KruegerTM` | Transverse Mercator forward/inverse using Karney's 2011 Krüger series — matches `GeographicLib::TransverseMercator` |
| `TransverseMercatorExact` | Exact Transverse Mercator via Jacobi elliptic functions — matches `GeographicLib::TransverseMercatorExact` |
| `KarneyUTMConverter` | Full UTM zone handling (forward/inverse) on top of both TM implementations |
| `EllipticFunction` | Jacobi elliptic functions and integrals, required by the exact TM |
| `DMS` | Degree/minute/second parsing and formatting — same as the JS `dms` module |

### Additional functionality — Snyder (independent)

The library also includes a classic UTM implementation based on **Snyder (1987)** that is fully independent of GeographicLib:

| Module | Source | Accuracy |
|---|---|---|
| `UTMConverter`, `UTMPoint`, `GeoPoint` | Snyder, J.P. (1987): *Map Projections — A Working Manual*, pp. 57–64. USGS Professional Paper 1395. | ~1 mm |

> Use `KarneyUTMConverter` when nanometer-level accuracy is required. Use `UTMConverter` for lightweight, dependency-free UTM conversions.

## Installation

Add to your `pubspec.yaml`:

```yaml
dependencies:
  geographiclib_dart:
    git:
      url: https://github.com/Fxztam/geographiclib-dart.git
```

Or, once published on pub.dev:

```yaml
dependencies:
  geographiclib_dart: ^2.0.0
```

## Usage

```dart
import 'package:geographiclib_dart/geographiclib.dart';

// Geodesic inverse: distance Berlin → Munich
final g = Geodesic.WGS84;
final r = g.inverse(52.52, 13.41, 48.14, 11.58);
print(r.s12); // ~504 103 m

// DMS decode / encode
final d = DMS.Decode("48d08'02.7\"N");
print(DMS.Encode(d.val, DMS.SECOND, 1, DMS.LATITUDE)); // 48°08'02.7"N

// UTM (Karney, ≤ 5 nm)
final utm = KarneyUTMConverter.fromLatLon(52.52, 13.41);
print('${utm.zone}${utm.hemisphere}  E=${utm.easting.toStringAsFixed(3)}');

// Transverse Mercator (Krüger series)
final tm = KruegerTM(); // WGS84 defaults
final fwd = tm.forward(9.0, 52.0, 13.41);  // central meridian, lat, lon
```

You can also import individual modules:

```dart
import 'package:geographiclib_dart/geodesic.dart';   // Geodesic only
import 'package:geographiclib_dart/dms.dart';          // DMS only
import 'package:geographiclib_dart/utm.dart';          // UTM only
```

## Running tests

```sh
git clone https://github.com/Fxztam/geographiclib-dart.git
cd geographiclib-dart
dart pub get
dart test
```

Expected output:

```
Building package executable... (1.9s)
Built test:test.
00:00 +271: All tests passed!
```

## Running the demo

A ready-to-run demo script is included in the `example/` folder:

```sh
git clone https://github.com/Fxztam/geographiclib-dart.git
cd geographiclib-dart/example
dart pub get
dart run bin/test-geolib.dart
```

Expected output:

```
=== Geodesic inverse: Berlin → Munich ===
Distance s12 = 504289.453 m
Azimuth  azi1 = -164.317446°

=== DMS decode / encode ===
Decoded value  = 48.134083333333336°
Re-encoded DMS = 48°08'02.7"N

=== UTM Karney (Berlin 52.52°N 13.41°E) ===
Zone:      33N
Easting:   392118.487 m
Northing:  5820064.675 m
γ (conv):  -1.261890431°
k (scale): 0.999742858842

=== KruegerTM forward (lon0=9°, lat=52°, lon=13.41°) ===
x (easting from CM):  302797.714552 m
y (northing):         5772534.263476 m
γ (convergence):      3.477748932°
k (scale):            1.001125279691
```

## Origins

Original C++/JS implementation by Charles Karney (MIT/X11 License):
- https://geographiclib.sourceforge.io/
- https://github.com/geographiclib/geographiclib-js

Dart port by Friedhold Matz (2026), MIT/X11 License.

## License

MIT — see [LICENSE](LICENSE).
