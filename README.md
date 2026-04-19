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
final utm = KarneyUTM.fromLatLon(52.52, 13.41);
print('${utm.zone}${utm.band}  E=${utm.easting.toStringAsFixed(3)}');

// Transverse Mercator (Krüger series)
final tm = KruegerTM.WGS84;
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
dart test
```

271 tests, all passing.

## Origins

Original C++/JS implementation by Charles Karney (MIT/X11 License):
- https://geographiclib.sourceforge.io/
- https://github.com/geographiclib/geographiclib-js

Dart port by Friedhold Matz (2026), MIT/X11 License.

## License

MIT — see [LICENSE](LICENSE).
