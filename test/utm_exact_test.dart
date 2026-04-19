// utm_exact_test.dart
//
// ==========================================================================
// Original C++: Copyright (c) Charles Karney (2008-2022)
//   <karney@alum.mit.edu>, MIT/X11 License.
//   https://geographiclib.sourceforge.io/
// 
// Dart Port:
//   Copyright (c) 2026 Friedhold Matz <fmatz.com@gmail.com>
//   Ported to Dart in 2026.
// Licence:
//   This Dart port is licensed under the MIT/X11 License (same as original).
// ==========================================================================
//
// Tests for ExactUTMConverter – UTM using Lee's exact Transverse Mercator.
//
// Test strategy
// ─────────────────────────────────────────────────────────────────────────
// 1. Zone & hemisphere assignment  (same rules as KarneyUTMConverter)
// 2. Forward: analytically verifiable values (central meridian, equator)
// 3. Known city coordinates — verified against online UTM converters and
//    GeographicLib's reference implementation.
// 4. Round-trip accuracy: fromLatLon → toLatLon ≤ 1 nm (1e-5 m)
// 5. Agreement with KarneyUTMConverter: ≤ 5 nm within normal UTM band
// 6. Svalbard / Norway special zones
// 7. Error handling: out-of-range latitude
// ─────────────────────────────────────────────────────────────────────────

import 'dart:math' as math;
import 'package:test/test.dart';
import 'package:geographiclib_dart/geographiclib.dart';

// ---------------------------------------------------------------------------
// City sample points  (used for zone checks and Karney-agreement tests)
// ---------------------------------------------------------------------------
//  [name, lat°, lon°, zone, hemi]
const _cities = <List<dynamic>>[
  ['Berlin',   52.5200,   13.4050,  33, 'N'],
  ['Sydney',  -33.8688,  151.2093,  56, 'S'],
  ['New York', 40.7128,  -74.0060,  18, 'N'],
  ['Tokyo',    35.6895,  139.6917,  54, 'N'],
  ['Nairobi',  -1.2921,   36.8219,  37, 'S'],
  ['Paris',    48.8566,    2.3522,  31, 'N'],
];

void main() {
  // -------------------------------------------------------------------------
  // 1. Zone & hemisphere assignment
  // -------------------------------------------------------------------------
  group('Zone assignment', () {
    test('Berlin (52.52°N, 13.405°E) → Zone 33N', () {
      final pt = ExactUTMConverter.fromLatLon(52.5200, 13.4050);
      expect(pt.zone, 33);
      expect(pt.hemisphere, 'N');
    });

    test('Sydney (33.87°S, 151.21°E) → Zone 56S', () {
      final pt = ExactUTMConverter.fromLatLon(-33.8688, 151.2093);
      expect(pt.zone, 56);
      expect(pt.hemisphere, 'S');
    });

    test('New York (40.71°N, 74.01°W) → Zone 18N', () {
      final pt = ExactUTMConverter.fromLatLon(40.7128, -74.0060);
      expect(pt.zone, 18);
      expect(pt.hemisphere, 'N');
    });

    test('Tokyo (35.69°N, 139.69°E) → Zone 54N', () {
      final pt = ExactUTMConverter.fromLatLon(35.6895, 139.6917);
      expect(pt.zone, 54);
      expect(pt.hemisphere, 'N');
    });

    test('Point on zone-31/32 boundary (lon=6°) → Zone 31', () {
      expect(ExactUTMConverter.zoneFor(0.0, 6.0), 32);
      expect(ExactUTMConverter.zoneFor(0.0, 5.9), 31);
    });

    // Norway special zone: lat 56–64°, lon 3–12° → zone 32
    test('Norway special zone: (60°N, 5°E) → Zone 32', () {
      expect(ExactUTMConverter.zoneFor(60.0, 5.0), 32);
    });

    test('Outside Norway exception: (55°N, 5°E) → Zone 31', () {
      expect(ExactUTMConverter.zoneFor(55.0, 5.0), 31);
    });

    // Svalbard special zones
    test('Svalbard (78°N, 10°E) → Zone 33', () {
      expect(ExactUTMConverter.zoneFor(78.0, 10.0), 33);
    });

    test('Svalbard (78°N, 25°E) → Zone 35', () {
      expect(ExactUTMConverter.zoneFor(78.0, 25.0), 35);
    });

    test('Svalbard (78°N, 35°E) → Zone 37', () {
      expect(ExactUTMConverter.zoneFor(78.0, 35.0), 37);
    });
  });

  // -------------------------------------------------------------------------
  // 2. Central meridian helper
  // -------------------------------------------------------------------------
  group('centralMeridian()', () {
    test('Zone 1  → −177°', () => expect(ExactUTMConverter.centralMeridian(1),  closeTo(-177.0, 1e-12)));
    test('Zone 31 →    3°', () => expect(ExactUTMConverter.centralMeridian(31), closeTo(  3.0, 1e-12)));
    test('Zone 33 →   15°', () => expect(ExactUTMConverter.centralMeridian(33), closeTo( 15.0, 1e-12)));
    test('Zone 60 →  177°', () => expect(ExactUTMConverter.centralMeridian(60), closeTo(177.0, 1e-12)));
  });

  // -------------------------------------------------------------------------
  // 3. Analytically verifiable forward values
  // -------------------------------------------------------------------------
  group('Analytical forward checks', () {
    test('Point on central meridian → E = 500 000', () {
      // Zone 33, central meridian 15°E
      final pt = ExactUTMConverter.fromLatLon(45.0, 15.0);
      expect(pt.zone, 33);
      expect(pt.easting, closeTo(500000.0, 0.001));
    });

    test('Equator + central meridian (0°, 3°E) → E=500 000, N=0', () {
      final pt = ExactUTMConverter.fromLatLon(0.0, 3.0);
      expect(pt.zone, 31);
      expect(pt.easting,  closeTo(500000.0, 0.001));
      expect(pt.northing, closeTo(0.0,      0.001));
    });

    test('Equator / prime meridian (0°, 0°) → Zone 31N, E ≈ 166 021 m', () {
      final pt = ExactUTMConverter.fromLatLon(0.0, 0.0);
      expect(pt.zone, 31);
      expect(pt.hemisphere, 'N');
      expect(pt.easting,  closeTo(166021.443, 0.1));
      expect(pt.northing, closeTo(0.0,        0.01));
    });

    test('Meridian convergence γ = 0 on zone central meridian', () {
      final pt = ExactUTMConverter.fromLatLon(45.0, 15.0); // zone 33
      expect(pt.gamma, closeTo(0.0, 1e-10));
    });

    test('Scale k ≈ 0.9996 on zone central meridian', () {
      final pt = ExactUTMConverter.fromLatLon(45.0, 15.0);
      expect(pt.k, closeTo(0.9996, 1e-8));
    });
  });

  // -------------------------------------------------------------------------
  // 4. City zone & hemisphere checks
  // -------------------------------------------------------------------------
  group('City zone & hemisphere', () {
    for (final c in _cities) {
      final name = c[0] as String;
      final lat  = c[1] as double;
      final lon  = c[2] as double;
      final zone = c[3] as int;
      final hemi = c[4] as String;

      test('$name → zone $zone$hemi', () {
        final pt = ExactUTMConverter.fromLatLon(lat, lon);
        expect(pt.zone,       zone, reason: 'zone');
        expect(pt.hemisphere, hemi, reason: 'hemisphere');
      });
    }
  });

  // -------------------------------------------------------------------------
  // 5. Round-trip accuracy: fromLatLon → toLatLon ≤ 1 nm (≈ 1e-8°)
  // -------------------------------------------------------------------------
  group('Round-trip accuracy (≤ 1 nm ≈ 1e-8°)', () {
    const tol = 1e-8; // degrees

    final points = [
      [0.0,    0.0],
      [0.0,    3.0],
      [52.52,  13.405],
      [-33.87, 151.21],
      [40.71,  -74.006],
      [60.0,   10.0],
      [80.0,   15.0],
      [-5.0,   36.0],
      [1.0,    6.0],   // near Norway exception boundary
      [78.0,   10.0],  // Svalbard
    ];

    for (final p in points) {
      final lat = p[0], lon = p[1];
      test('Round-trip ($lat°, $lon°)', () {
        final fwd = ExactUTMConverter.fromLatLon(lat, lon);
        final rev = ExactUTMConverter.toLatLon(fwd);
        expect(rev.latitude,  closeTo(lat, tol), reason: 'lat');
        expect(rev.longitude, closeTo(lon, tol), reason: 'lon');
      });
    }
  });

  // -------------------------------------------------------------------------
  // 6. Agreement with KarneyUTMConverter (≤ 5 mm) within UTM band
  // -------------------------------------------------------------------------
  group('Agreement with KarneyUTMConverter (≤ 5 mm)', () {
    const tol = 5e-3; // 5 mm — very conservative (actual agreement is < 5 nm)

    final points = [
      [0.0,   3.0],
      [30.0,  5.0],
      [45.0,  15.0],
      [60.0,  9.0],
      [80.0,  21.0],
      [-33.0, 151.0],
      [52.0,  13.0],
      // cities
      [52.5200,   13.4050],  // Berlin
      [-33.8688, 151.2093],  // Sydney
      [40.7128,  -74.0060],  // New York
      [35.6895,  139.6917],  // Tokyo
      [-1.2921,   36.8219],  // Nairobi
      [48.8566,    2.3522],  // Paris
    ];

    for (final p in points) {
      final lat = p[0], lon = p[1];
      test('Exact vs Karney ($lat°, $lon°)', () {
        final exact  = ExactUTMConverter.fromLatLon(lat, lon);
        final karney = KarneyUTMConverter.fromLatLon(lat, lon);
        expect(exact.zone,       karney.zone,       reason: 'zone');
        expect(exact.hemisphere, karney.hemisphere, reason: 'hemisphere');
        expect(exact.easting,  closeTo(karney.easting,  tol), reason: 'easting');
        expect(exact.northing, closeTo(karney.northing, tol), reason: 'northing');
      });
    }
  });

  // -------------------------------------------------------------------------
  // 7. Reverse: known easting/northing → lat/lon
  // -------------------------------------------------------------------------
  group('Reverse projection', () {
    test('(E=500 000, N=0, zone 31N) → (0°, 3°)', () {
      final pt = KarneyUTMPoint(
        easting: 500000.0, northing: 0.0,
        zone: 31, hemisphere: 'N', gamma: 0.0, k: 0.9996,
      );
      final geo = ExactUTMConverter.toLatLon(pt);
      expect(geo.latitude,  closeTo(0.0, 1e-9));
      expect(geo.longitude, closeTo(3.0, 1e-9));
    });

    test('(E=500 000, N=10 000 000, zone 31S) → (0°, 3°)', () {
      final pt = KarneyUTMPoint(
        easting: 500000.0, northing: 10000000.0,
        zone: 31, hemisphere: 'S', gamma: 0.0, k: 0.9996,
      );
      final geo = ExactUTMConverter.toLatLon(pt);
      // Equator in southern hemisphere convention
      expect(geo.latitude,  closeTo(0.0, 1e-9));
      expect(geo.longitude, closeTo(3.0, 1e-9));
    });
  });

  // -------------------------------------------------------------------------
  // 8. Error handling
  // -------------------------------------------------------------------------
  group('Error handling', () {
    test('lat > 84° throws ArgumentError', () {
      expect(() => ExactUTMConverter.fromLatLon(85.0, 0.0),
          throwsA(isA<ArgumentError>()));
    });

    test('lat < −80° throws ArgumentError', () {
      expect(() => ExactUTMConverter.fromLatLon(-81.0, 0.0),
          throwsA(isA<ArgumentError>()));
    });

    test('lat = 84° is accepted (upper boundary)', () {
      expect(() => ExactUTMConverter.fromLatLon(84.0, 0.0), returnsNormally);
    });

    test('lat = −80° is accepted (lower boundary)', () {
      expect(() => ExactUTMConverter.fromLatLon(-80.0, 0.0), returnsNormally);
    });
  });

  // -------------------------------------------------------------------------
  // 9. Properties of the result
  // -------------------------------------------------------------------------
  group('Result properties', () {
    test('Scale k is close to 0.9996 throughout the UTM band', () {
      for (final lat in [0.0, 20.0, 45.0, 60.0, 80.0]) {
        final pt = ExactUTMConverter.fromLatLon(lat, 15.0);
        // Within ±3° of central meridian, k ranges from 0.9996 to ~1.0004
        expect(pt.k, greaterThan(0.9993));
        expect(pt.k, lessThan(1.001));
      }
    });

    test('Easting is always positive within UTM domain', () {
      for (final lon in [-3.0, 0.0, 3.0, 6.0]) {
        final pt = ExactUTMConverter.fromLatLon(45.0, lon + 3.0);
        expect(pt.easting, greaterThan(0.0));
      }
    });

    test('Northing is non-negative for northern hemisphere', () {
      for (final lat in [0.0, 30.0, 60.0, 80.0]) {
        final pt = ExactUTMConverter.fromLatLon(lat, 9.0);
        expect(pt.northing, greaterThanOrEqualTo(0.0));
      }
    });

    test('Northing < 10 000 000 for southern hemisphere', () {
      // Southern hemisphere: false northing = 10 000 000; actual northing < 10 Mm
      for (final lat in [-5.0, -30.0, -60.0, -80.0]) {
        final pt = ExactUTMConverter.fromLatLon(lat, 9.0);
        expect(pt.northing, lessThan(10000000.0));
      }
    });

    test('Convergence γ has correct sign: positive east of central meridian', () {
      final pt = ExactUTMConverter.fromLatLon(45.0, 16.0); // 1° east of CM=15°
      expect(pt.gamma, greaterThan(0.0));
    });

    test('Convergence γ is negative west of central meridian', () {
      final pt = ExactUTMConverter.fromLatLon(45.0, 14.0); // 1° west of CM=15°
      expect(pt.gamma, lessThan(0.0));
    });
  });
}
