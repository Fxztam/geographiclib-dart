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

import 'package:test/test.dart';
import 'package:geographiclib_dart/geographiclib.dart';

void main() {
  const double tol = 0.001; // 1 mm tolerance for coordinates

  // -------------------------------------------------------------------------
  // UTM forward (LatLon → UTM) – zone / hemisphere assignment
  // -------------------------------------------------------------------------
  group('UTMConverter.fromLatLon – zone & hemisphere', () {
    test('Berlin (52.5200°N, 13.4050°E) → Zone 33N', () {
      final utm = UTMConverter.fromLatLon(52.5200, 13.4050);
      expect(utm.zone, equals(33));
      expect(utm.hemisphere, equals('N'));
    });

    test('Sydney (−33.8688°, 151.2093°E) → Zone 56S', () {
      final utm = UTMConverter.fromLatLon(-33.8688, 151.2093);
      expect(utm.zone, equals(56));
      expect(utm.hemisphere, equals('S'));
    });

    test('New York (40.7128°N, −74.0060°W) → Zone 18N', () {
      final utm = UTMConverter.fromLatLon(40.7128, -74.0060);
      expect(utm.zone, equals(18));
      expect(utm.hemisphere, equals('N'));
    });

    test('Tokyo (35.6895°N, 139.6917°E) → Zone 54N', () {
      final utm = UTMConverter.fromLatLon(35.6895, 139.6917);
      expect(utm.zone, equals(54));
      expect(utm.hemisphere, equals('N'));
    });

    test('Latitude out of range throws ArgumentError', () {
      expect(() => UTMConverter.fromLatLon(85.0, 0.0),
          throwsA(isA<ArgumentError>()));
      expect(() => UTMConverter.fromLatLon(-81.0, 0.0),
          throwsA(isA<ArgumentError>()));
    });
  });

  // -------------------------------------------------------------------------
  // UTM forward – analytically verifiable values
  // -------------------------------------------------------------------------
  group('UTMConverter.fromLatLon – known values', () {
    // A point on the central meridian of a zone always has E = 500 000 m.
    test('Point on central meridian → E = 500 000', () {
      // Zone 33, central meridian = 15°E
      final utm = UTMConverter.fromLatLon(45.0, 15.0);
      expect(utm.zone, equals(33));
      expect(utm.easting, closeTo(500000.0, tol));
    });

    test('Equator / prime meridian (0°, 0°) → Zone 31N, E ≈ 166 021', () {
      final utm = UTMConverter.fromLatLon(0.0, 0.0);
      expect(utm.zone, equals(31));
      expect(utm.hemisphere, equals('N'));
      // 3° west of the central meridian (3°E): roughly 500000 - 333979 = 166021
      expect(utm.easting, closeTo(166021.443, tol));
      expect(utm.northing, closeTo(0.0, tol));
    });

    test('Equator / central meridian (0°, 3°E) → E = 500 000, N = 0', () {
      final utm = UTMConverter.fromLatLon(0.0, 3.0);
      expect(utm.zone, equals(31));
      expect(utm.easting, closeTo(500000.0, tol));
      expect(utm.northing, closeTo(0.0, tol));
    });
  });

  // -------------------------------------------------------------------------
  // UTM inverse (UTM → LatLon)
  // -------------------------------------------------------------------------
  group('UTMConverter.toLatLon', () {
    test('Zone 33N Berlin round-trip', () {
      final src = GeoPoint(52.5200, 13.4050);
      final utm = UTMConverter.fromLatLon(src.latitude, src.longitude);
      final geo = UTMConverter.toLatLon(utm);
      expect(geo.latitude, closeTo(src.latitude, 1e-7));
      expect(geo.longitude, closeTo(src.longitude, 1e-7));
    });

    test('Zone 56S Sydney round-trip', () {
      final src = GeoPoint(-33.8688, 151.2093);
      final utm = UTMConverter.fromLatLon(src.latitude, src.longitude);
      final geo = UTMConverter.toLatLon(utm);
      expect(geo.latitude, closeTo(src.latitude, 1e-7));
      expect(geo.longitude, closeTo(src.longitude, 1e-7));
    });

    test('Zone 18N New York round-trip', () {
      final src = GeoPoint(40.7128, -74.0060);
      final utm = UTMConverter.fromLatLon(src.latitude, src.longitude);
      final geo = UTMConverter.toLatLon(utm);
      expect(geo.latitude, closeTo(src.latitude, 1e-7));
      expect(geo.longitude, closeTo(src.longitude, 1e-7));
    });
  });

  // -------------------------------------------------------------------------
  // DMS ↔ decimal
  // -------------------------------------------------------------------------
  group('DMS.fromDecimal', () {
    test('48.8566 → 48°51\'23.76"', () {
      final d = DMS.fromDecimal(48.8566);
      expect(d.degrees, equals(48));
      expect(d.minutes, equals(51));
      expect(d.seconds, closeTo(23.76, 0.01));
    });

    test('negative value −33.8688 keeps sign on degrees', () {
      final d = DMS.fromDecimal(-33.8688);
      expect(d.degrees, equals(-33));
      expect(d.minutes, equals(52));
      expect(d.seconds, closeTo(7.68, 0.01));
    });
  });

  group('DMS.toDecimal', () {
    test('48 51 23.76 N → 48.8566...', () {
      final dec = DMS.toDecimal(48, 51, 23.76, 'N');
      expect(dec, closeTo(48.8566, 0.0001));
    });

    test('13 24 18.00 E → 13.405', () {
      final dec = DMS.toDecimal(13, 24, 18.0, 'E');
      expect(dec, closeTo(13.405, 0.0001));
    });

    test('S hemisphere returns negative', () {
      final dec = DMS.toDecimal(33, 52, 7.68, 'S');
      expect(dec, closeTo(-33.8688, 0.0001));
    });

    test('W hemisphere returns negative', () {
      final dec = DMS.toDecimal(74, 0, 21.6, 'W');
      expect(dec, closeTo(-74.006, 0.0001));
    });

    test('round-trip: decimal → DMS → decimal', () {
      const double src = 52.5200;
      final dms = DMS.fromDecimal(src);
      final dec = DMS.toDecimal(dms.degrees, dms.minutes, dms.seconds, 'N');
      expect(dec, closeTo(src, 1e-9));
    });
  });

  group('DMS.format', () {
    test('latitude 48.8566 → "48°51\'23.76\"N"', () {
      expect(DMS.format(48.8566), equals("48°51'23.76\"N"));
    });

    test('longitude 13.405', () {
      expect(DMS.format(13.405, isLatitude: false), equals("13°24'18.00\"E"));
    });

    test('negative longitude', () {
      expect(DMS.format(-74.006, isLatitude: false), equals("74°0'21.60\"W"));
    });
  });

  group('DMS.parse', () {
    test('standard DMS string', () {
      expect(DMS.parse("48°51'23.76\"N"), closeTo(48.8566, 0.0001));
    });

    test('colon-separated string', () {
      expect(DMS.parse('48:51:23.76N'), closeTo(48.8566, 0.0001));
    });

    test('space-separated string', () {
      expect(DMS.parse('48 51 23.76 N'), closeTo(48.8566, 0.0001));
    });

    test('south hemisphere', () {
      expect(DMS.parse("33°52'7.68\"S"), closeTo(-33.8688, 0.0001));
    });

    test('invalid string throws FormatException', () {
      expect(() => DMS.parse('not a coordinate'), throwsA(isA<FormatException>()));
    });
  });
}
