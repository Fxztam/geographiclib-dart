# geographiclib_dart

Dart port of [Charles Karney's GeographicLib](https://geographiclib.sourceforge.io/).

A pure-Dart library — works in any Dart project and in Flutter apps.

## Features

| Class / API | Description | Accuracy |
|---|---|---|
| `Geodesic`, `GeodesicLine`, `PolygonArea` | Geodesic direct, inverse, area on WGS84 ellipsoid | < 15 nm |
| `DMS` | DMS string parsing & formatting (Karney API + legacy helpers) | — |
| `KruegerTM`, `TransverseMercatorResult` | Transverse Mercator — 6th-order Krüger series | ≤ 5 nm |
| `TransverseMercatorExact` | Transverse Mercator — Lee's exact method | exact |
| `EllipticFunction` | Jacobi elliptic functions and integrals | — |
| `KarneyUTM`, `UTMResult` | UTM ↔ geographic (Karney) | ≤ 5 nm |
| `UTMConverter`, `UTMPoint`, `GeoPoint` | UTM ↔ geographic (Snyder) | ~1 mm |

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
  geographiclib_dart: ^1.0.0
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
