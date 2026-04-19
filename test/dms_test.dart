// dms_test.dart
// Dart port of dmstest.js from geographiclib-js by Charles Karney.
// https://github.com/geographiclib/geographiclib-js/blob/main/dms/test/dmstest.js
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

import 'dart:math' as math;
import 'package:test/test.dart';
import 'package:geographiclib_dart/geographiclib.dart';

// Helpers
bool _decodedEqual(({double val, int ind}) a, ({double val, int ind}) b) {
  if (a.ind != b.ind) return false;
  // Handle NaN == NaN comparison
  if (a.val.isNaN && b.val.isNaN) return true;
  return a.val == b.val;
}

void main() {
  group('DMSTest', () {
    // -----------------------------------------------------------------------
    // Decode
    // -----------------------------------------------------------------------
    test('check decode', () {
      expect(
        _decodedEqual(DMS.Decode('E7:33:36'), DMS.Decode('-7.56W')),
        isTrue,
      );
    });

    // -----------------------------------------------------------------------
    // Encode
    // -----------------------------------------------------------------------
    test('check encode — basic', () {
      expect(DMS.Encode(-7.56, DMS.DEGREE, 2), equals('-7.56°'));
      expect(DMS.Encode(-7.56, DMS.MINUTE, 1), equals("-7°33.6'"));
      expect(DMS.Encode(-7.56, DMS.SECOND, 0), equals('-7°33\'36"'));

      expect(DMS.Encode(-7.56, DMS.DEGREE, 2, DMS.NONE, ':'), equals('-7.56'));
      expect(DMS.Encode(-7.56, DMS.MINUTE, 1, DMS.NONE, ':'), equals('-7:33.6'));
      expect(DMS.Encode(-7.56, DMS.SECOND, 0, DMS.NONE, ':'), equals('-7:33:36'));

      expect(DMS.Encode(-7.56, DMS.DEGREE, 2, DMS.LATITUDE), equals('07.56°S'));
      expect(DMS.Encode(-7.56, DMS.MINUTE, 1, DMS.LATITUDE), equals("07°33.6'S"));
      expect(DMS.Encode(-7.56, DMS.SECOND, 0, DMS.LATITUDE), equals('07°33\'36"S'));

      expect(DMS.Encode(-7.56, DMS.DEGREE, 2, DMS.LATITUDE, ':'), equals('07.56S'));
      expect(DMS.Encode(-7.56, DMS.MINUTE, 1, DMS.LATITUDE, ':'), equals('07:33.6S'));
      expect(DMS.Encode(-7.56, DMS.SECOND, 0, DMS.LATITUDE, ':'), equals('07:33:36S'));
    });

    test('check encode — zero-fill', () {
      // t = -(1 + 2/60 + 2.99/3600)  ≈  -1°02'02.99"
      final double t = -(1 + 2 / 60 + 2.99 / 3600);

      expect(DMS.Encode( t, DMS.DEGREE, 0, DMS.NONE     ), equals('-1°'));
      expect(DMS.Encode( t, DMS.DEGREE, 0, DMS.LATITUDE ), equals('01°S'));
      expect(DMS.Encode( t, DMS.DEGREE, 0, DMS.LONGITUDE), equals('001°W'));
      expect(DMS.Encode(-t, DMS.DEGREE, 0, DMS.AZIMUTH  ), equals('001°'));

      expect(DMS.Encode( t, DMS.DEGREE, 1, DMS.NONE     ), equals('-1.0°'));
      expect(DMS.Encode( t, DMS.DEGREE, 1, DMS.LATITUDE ), equals('01.0°S'));
      expect(DMS.Encode( t, DMS.DEGREE, 1, DMS.LONGITUDE), equals('001.0°W'));
      expect(DMS.Encode(-t, DMS.DEGREE, 1, DMS.AZIMUTH  ), equals('001.0°'));

      expect(DMS.Encode( t, DMS.MINUTE, 0, DMS.NONE     ), equals("-1°02'"));
      expect(DMS.Encode( t, DMS.MINUTE, 0, DMS.LATITUDE ), equals("01°02'S"));
      expect(DMS.Encode( t, DMS.MINUTE, 0, DMS.LONGITUDE), equals("001°02'W"));
      expect(DMS.Encode(-t, DMS.MINUTE, 0, DMS.AZIMUTH  ), equals("001°02'"));

      expect(DMS.Encode( t, DMS.MINUTE, 1, DMS.NONE     ), equals("-1°02.0'"));
      expect(DMS.Encode( t, DMS.MINUTE, 1, DMS.LATITUDE ), equals("01°02.0'S"));
      expect(DMS.Encode( t, DMS.MINUTE, 1, DMS.LONGITUDE), equals("001°02.0'W"));
      expect(DMS.Encode(-t, DMS.MINUTE, 1, DMS.AZIMUTH  ), equals("001°02.0'"));

      expect(DMS.Encode( t, DMS.SECOND, 0, DMS.NONE     ), equals('-1°02\'03"'));
      expect(DMS.Encode( t, DMS.SECOND, 0, DMS.LATITUDE ), equals('01°02\'03"S'));
      expect(DMS.Encode( t, DMS.SECOND, 0, DMS.LONGITUDE), equals('001°02\'03"W'));
      expect(DMS.Encode(-t, DMS.SECOND, 0, DMS.AZIMUTH  ), equals('001°02\'03"'));

      expect(DMS.Encode( t, DMS.SECOND, 1, DMS.NONE     ), equals('-1°02\'03.0"'));
      expect(DMS.Encode( t, DMS.SECOND, 1, DMS.LATITUDE ), equals('01°02\'03.0"S'));
      expect(DMS.Encode( t, DMS.SECOND, 1, DMS.LONGITUDE), equals('001°02\'03.0"W'));
      expect(DMS.Encode(-t, DMS.SECOND, 1, DMS.AZIMUTH  ), equals('001°02\'03.0"'));
    });

    // -----------------------------------------------------------------------
    // Decode special values
    // -----------------------------------------------------------------------
    test('check decode special', () {
      expect(DMS.Decode(' +0 ').val,  equals(0.0));
      // -0 must decode as negative zero
      expect(1 / DMS.Decode('-0  ').val, equals(double.negativeInfinity));
      expect(DMS.Decode(' nan').val,  isNaN);
      expect(DMS.Decode('+inf').val,  equals(double.infinity));
      expect(DMS.Decode(' inf').val,  equals(double.infinity));
      expect(DMS.Decode('-inf').val,  equals(double.negativeInfinity));

      expect(DMS.Decode(' +0N').val,  equals(0.0));
      // -0N  → val = -0
      expect(1 / DMS.Decode('-0N ').val, equals(double.negativeInfinity));
      // +0S  → sign flipped → val = -0
      expect(1 / DMS.Decode('+0S ').val, equals(double.negativeInfinity));
      // -0S  → double negation → val = +0
      expect(1 / DMS.Decode(' -0S').val, equals(double.infinity));
    });

    // -----------------------------------------------------------------------
    // Encode rounding
    // -----------------------------------------------------------------------
    test('check encode rounding', () {
      // Dart's toStringAsFixed rounds ties away from zero (like JS .toFixed)
      expect(DMS.Encode(double.nan,          DMS.DEGREE, 0), equals('nan'));
      expect(DMS.Encode(double.infinity,     DMS.DEGREE, 0), equals('inf'));
      expect(DMS.Encode(double.negativeInfinity, DMS.DEGREE, 0), equals('-inf'));

      expect(DMS.Encode(-3.5, DMS.DEGREE, 0), equals('-4°'));
      expect(DMS.Encode(-2.5, DMS.DEGREE, 0), equals('-3°')); // round-half-up
      expect(DMS.Encode(-1.5, DMS.DEGREE, 0), equals('-2°'));
      expect(DMS.Encode(-0.5, DMS.DEGREE, 0), equals('-1°')); // round-half-up
      expect(DMS.Encode(-0.0, DMS.DEGREE, 0), equals('-0°'));
      expect(DMS.Encode( 0.0, DMS.DEGREE, 0), equals( '0°'));
      expect(DMS.Encode( 0.5, DMS.DEGREE, 0), equals( '1°')); // round-half-up
      expect(DMS.Encode( 1.5, DMS.DEGREE, 0), equals( '2°'));
      expect(DMS.Encode( 2.5, DMS.DEGREE, 0), equals( '3°')); // round-half-up
      expect(DMS.Encode( 3.5, DMS.DEGREE, 0), equals( '4°'));

      expect(DMS.Encode(-1.75, DMS.DEGREE, 1), equals('-1.8°'));
      expect(DMS.Encode(-1.25, DMS.DEGREE, 1), equals('-1.3°')); // round-half-up
      expect(DMS.Encode(-0.75, DMS.DEGREE, 1), equals('-0.8°'));
      expect(DMS.Encode(-0.25, DMS.DEGREE, 1), equals('-0.3°')); // round-half-up
      expect(DMS.Encode(-0.0,  DMS.DEGREE, 1), equals('-0.0°'));
      expect(DMS.Encode( 0.0,  DMS.DEGREE, 1), equals( '0.0°'));
      expect(DMS.Encode( 0.25, DMS.DEGREE, 1), equals( '0.3°')); // round-half-up
      expect(DMS.Encode( 0.75, DMS.DEGREE, 1), equals( '0.8°'));
      expect(DMS.Encode( 1.25, DMS.DEGREE, 1), equals( '1.3°')); // round-half-up
      expect(DMS.Encode( 1.75, DMS.DEGREE, 1), equals( '1.8°'));

      expect(DMS.Encode(1e20, DMS.DEGREE, 0),
          equals('100000000000000000000°'));
      expect(DMS.Encode(1e21, DMS.DEGREE, 0), equals('1e21'));
    });
  });
}
