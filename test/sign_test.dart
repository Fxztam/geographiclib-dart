// sign_test.dart
// Dart port of signtest.js from geographiclib-geodesic v2.2.
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
// Tests signed-zero and special-value behaviour of GeoMath functions.

import 'package:test/test.dart';
import 'package:geographiclib_dart/geodesic.dart';

// ---------------------------------------------------------------------------
// Helper: 'Object.is' semantics for double
// (distinguishes +0.0 / -0.0 and treats NaN == NaN)
// ---------------------------------------------------------------------------
void expectExact(double actual, double expected, [String? reason]) {
  if (expected.isNaN) {
    expect(actual.isNaN, isTrue,
        reason: reason ?? 'Expected NaN but got $actual');
  } else if (expected == 0.0) {
    expect(actual == 0.0, isTrue,
        reason: reason ?? 'Expected ±0.0 but got $actual');
    expect(actual.isNegative, equals(expected.isNegative),
        reason: reason ??
            'Wrong sign of zero: expected '
                '${expected.isNegative ? "-0.0" : "+0.0"} but got $actual');
  } else {
    expect(actual, equals(expected), reason: reason);
  }
}

// Helper for sincosd: checks both s and c with exact semantics.
void expectSinCos(
    (double s, double c) v, double expectedS, double expectedC) {
  expectExact(v.$1, expectedS, 'sincosd s');
  expectExact(v.$2, expectedC, 'sincosd c');
}

void main() {
  const double eps = GeoMath.epsilon; // 2^-52

  group('SignTest', () {
    test('check AngRound', () {
      expectExact(GeoMath.angRound(-eps / 32), -eps / 32);
      expectExact(GeoMath.angRound(-eps / 64), -0.0);
      expectExact(GeoMath.angRound(-0.0), -0.0);
      expectExact(GeoMath.angRound(0.0), 0.0);
      expectExact(GeoMath.angRound(eps / 64), 0.0);
      expectExact(GeoMath.angRound(eps / 32), eps / 32);
      expectExact(GeoMath.angRound((1 - 2 * eps) / 64), (1 - 2 * eps) / 64);
      expectExact(GeoMath.angRound((1 - eps) / 64), 1 / 64);
      expectExact(GeoMath.angRound((1 - eps / 2) / 64), 1 / 64);
      expectExact(GeoMath.angRound((1 - eps / 4) / 64), 1 / 64);
      expectExact(GeoMath.angRound(1 / 64), 1 / 64);
      expectExact(GeoMath.angRound((1 + eps / 2) / 64), 1 / 64);
      expectExact(GeoMath.angRound((1 + eps) / 64), 1 / 64);
      expectExact(GeoMath.angRound((1 + 2 * eps) / 64), (1 + 2 * eps) / 64);
      expectExact(GeoMath.angRound((1 - eps) / 32), (1 - eps) / 32);
      expectExact(GeoMath.angRound((1 - eps / 2) / 32), 1 / 32);
      expectExact(GeoMath.angRound((1 - eps / 4) / 32), 1 / 32);
      expectExact(GeoMath.angRound(1 / 32), 1 / 32);
      expectExact(GeoMath.angRound((1 + eps / 2) / 32), 1 / 32);
      expectExact(GeoMath.angRound((1 + eps) / 32), (1 + eps) / 32);
      expectExact(GeoMath.angRound((1 - eps) / 16), (1 - eps) / 16);
      expectExact(GeoMath.angRound((1 - eps / 2) / 16), (1 - eps / 2) / 16);
      expectExact(GeoMath.angRound((1 - eps / 4) / 16), 1 / 16);
      expectExact(GeoMath.angRound(1 / 16), 1 / 16);
      expectExact(GeoMath.angRound((1 + eps / 4) / 16), 1 / 16);
      expectExact(GeoMath.angRound((1 + eps / 2) / 16), 1 / 16);
      expectExact(GeoMath.angRound((1 + eps) / 16), (1 + eps) / 16);
      expectExact(GeoMath.angRound((1 - eps) / 8), (1 - eps) / 8);
      expectExact(GeoMath.angRound((1 - eps / 2) / 8), (1 - eps / 2) / 8);
      expectExact(GeoMath.angRound((1 - eps / 4) / 8), 1 / 8);
      expectExact(GeoMath.angRound((1 + eps / 2) / 8), 1 / 8);
      expectExact(GeoMath.angRound((1 + eps) / 8), (1 + eps) / 8);
      expectExact(GeoMath.angRound(1 - eps), 1 - eps);
      expectExact(GeoMath.angRound(1 - eps / 2), 1 - eps / 2);
      expectExact(GeoMath.angRound(1 - eps / 4), 1.0);
      expectExact(GeoMath.angRound(1.0), 1.0);
      expectExact(GeoMath.angRound(1 + eps / 4), 1.0);
      expectExact(GeoMath.angRound(1 + eps / 2), 1.0);
      expectExact(GeoMath.angRound(1 + eps), 1 + eps);
      expectExact(GeoMath.angRound(90 - 64 * eps), 90 - 64 * eps);
      expectExact(GeoMath.angRound(90 - 32 * eps), 90.0);
      expectExact(GeoMath.angRound(90.0), 90.0);
    });

    test('check sincosd', () {
      const nan = double.nan;
      const inf = double.infinity;

      expectSinCos(GeoMath.sincosd(-inf), nan, nan);
      expectSinCos(GeoMath.sincosd(-810), -1.0, 0.0);
      expectSinCos(GeoMath.sincosd(-720), -0.0, 1.0);
      expectSinCos(GeoMath.sincosd(-630), 1.0, 0.0);
      expectSinCos(GeoMath.sincosd(-540), -0.0, -1.0);
      expectSinCos(GeoMath.sincosd(-450), -1.0, 0.0);
      expectSinCos(GeoMath.sincosd(-360), -0.0, 1.0);
      expectSinCos(GeoMath.sincosd(-270), 1.0, 0.0);
      expectSinCos(GeoMath.sincosd(-180), -0.0, -1.0);
      expectSinCos(GeoMath.sincosd(-90), -1.0, 0.0);
      expectSinCos(GeoMath.sincosd(-0.0), -0.0, 1.0);
      expectSinCos(GeoMath.sincosd(0.0), 0.0, 1.0);
      expectSinCos(GeoMath.sincosd(90), 1.0, 0.0);
      expectSinCos(GeoMath.sincosd(180), 0.0, -1.0);
      expectSinCos(GeoMath.sincosd(270), -1.0, 0.0);
      expectSinCos(GeoMath.sincosd(360), 0.0, 1.0);
      expectSinCos(GeoMath.sincosd(450), 1.0, 0.0);
      expectSinCos(GeoMath.sincosd(540), 0.0, -1.0);
      expectSinCos(GeoMath.sincosd(630), -1.0, 0.0);
      expectSinCos(GeoMath.sincosd(720), 0.0, 1.0);
      expectSinCos(GeoMath.sincosd(810), 1.0, 0.0);
      expectSinCos(GeoMath.sincosd(inf), nan, nan);
      expectSinCos(GeoMath.sincosd(nan), nan, nan);

      final v1 = GeoMath.sincosd(9.0);
      final v2 = GeoMath.sincosd(81.0);
      final v3 = GeoMath.sincosd(-123456789.0);
      expect(v1.$1, equals(v2.$2), reason: 'sin(9)==cos(81)');
      expect(v1.$1, equals(v3.$1), reason: 'sin(9)==sin(-123456789)');
      expect(v1.$2, equals(v2.$1), reason: 'cos(9)==sin(81)');
      expect(v1.$2, equals(-v3.$2), reason: 'cos(9)==-cos(-123456789)');
    });

    test('check atan2d', () {
      const nan = double.nan;
      const inf = double.infinity;

      expectExact(GeoMath.atan2d(0.0, -0.0), 180.0);
      expectExact(GeoMath.atan2d(-0.0, -0.0), -180.0);
      expectExact(GeoMath.atan2d(0.0, 0.0), 0.0);
      expectExact(GeoMath.atan2d(-0.0, 0.0), -0.0);
      expectExact(GeoMath.atan2d(0.0, -1.0), 180.0);
      expectExact(GeoMath.atan2d(-0.0, -1.0), -180.0);
      expectExact(GeoMath.atan2d(0.0, 1.0), 0.0);
      expectExact(GeoMath.atan2d(-0.0, 1.0), -0.0);
      expectExact(GeoMath.atan2d(-1.0, 0.0), -90.0);
      expectExact(GeoMath.atan2d(-1.0, -0.0), -90.0);
      expectExact(GeoMath.atan2d(1.0, 0.0), 90.0);
      expectExact(GeoMath.atan2d(1.0, -0.0), 90.0);
      expectExact(GeoMath.atan2d(1.0, -inf), 180.0);
      expectExact(GeoMath.atan2d(-1.0, -inf), -180.0);
      expectExact(GeoMath.atan2d(1.0, inf), 0.0);
      expectExact(GeoMath.atan2d(-1.0, inf), -0.0);
      expectExact(GeoMath.atan2d(inf, 1.0), 90.0);
      expectExact(GeoMath.atan2d(inf, -1.0), 90.0);
      expectExact(GeoMath.atan2d(-inf, 1.0), -90.0);
      expectExact(GeoMath.atan2d(-inf, -1.0), -90.0);
      expectExact(GeoMath.atan2d(inf, -inf), 135.0);
      expectExact(GeoMath.atan2d(-inf, -inf), -135.0);
      expectExact(GeoMath.atan2d(inf, inf), 45.0);
      expectExact(GeoMath.atan2d(-inf, inf), -45.0);
      expectExact(GeoMath.atan2d(nan, 1.0), nan);
      expectExact(GeoMath.atan2d(1.0, nan), nan);
      // atan2d(NaN, -1) == 180 - atan2d(NaN, 1)  (both NaN, Object.is)
      expectExact(GeoMath.atan2d(nan, -1.0), 180.0 - GeoMath.atan2d(nan, 1.0));
    });

    test('check sum', () {
      expectExact(GeoMath.sum(9.0, -9.0).$1, 0.0);
      expectExact(GeoMath.sum(-9.0, 9.0).$1, 0.0);
      expectExact(GeoMath.sum(-0.0, 0.0).$1, 0.0);
      expectExact(GeoMath.sum(0.0, -0.0).$1, 0.0);
      expectExact(GeoMath.sum(-0.0, -0.0).$1, -0.0);
      expectExact(GeoMath.sum(0.0, 0.0).$1, 0.0);
    });

    test('check AngNormalize', () {
      expectExact(GeoMath.angNormalize(-900.0), -180.0);
      expectExact(GeoMath.angNormalize(-720.0), -0.0);
      expectExact(GeoMath.angNormalize(-540.0), -180.0);
      expectExact(GeoMath.angNormalize(-360.0), -0.0);
      expectExact(GeoMath.angNormalize(-180.0), -180.0);
      expectExact(GeoMath.angNormalize(-0.0), -0.0);
      expectExact(GeoMath.angNormalize(0.0), 0.0);
      expectExact(GeoMath.angNormalize(180.0), 180.0);
      expectExact(GeoMath.angNormalize(360.0), 0.0);
      expectExact(GeoMath.angNormalize(540.0), 180.0);
      expectExact(GeoMath.angNormalize(720.0), 0.0);
      expectExact(GeoMath.angNormalize(900.0), 180.0);
    });

    test('check AngDiff', () {
      expectExact(GeoMath.angDiff(0.0, 0.0).$1, 0.0);
      expectExact(GeoMath.angDiff(0.0, -0.0).$1, -0.0);
      expectExact(GeoMath.angDiff(-0.0, 0.0).$1, 0.0);
      expectExact(GeoMath.angDiff(-0.0, -0.0).$1, 0.0);
      expectExact(GeoMath.angDiff(5.0, 365.0).$1, 0.0);
      expectExact(GeoMath.angDiff(365.0, 5.0).$1, -0.0);
      expectExact(GeoMath.angDiff(5.0, 185.0).$1, 180.0);
      expectExact(GeoMath.angDiff(185.0, 5.0).$1, -180.0);
      expectExact(GeoMath.angDiff(eps, 180.0).$1, 180.0);
      expectExact(GeoMath.angDiff(-eps, 180.0).$1, -180.0);
      expectExact(GeoMath.angDiff(eps, -180.0).$1, 180.0);
      expectExact(GeoMath.angDiff(-eps, -180.0).$1, -180.0);
      final x = 138.0 + 128.0 * eps;
      const y = -164.0;
      expectExact(GeoMath.angDiff(x, y).$1, 58.0 - 128.0 * eps);
    });

    test('azimuth with coincident point on equator', () {
      final geod = Geodesic.WGS84;
      // [lat1, lat2, expectedAzi]
      final cases = [
        [0.0, -0.0, 180.0],
        [-0.0, 0.0, 0.0],
      ];
      for (final c in cases) {
        final inv = geod.inverse(c[0], 0.0, c[1], 0.0);
        expectExact(inv.azi1, c[2]);
        expectExact(inv.azi2, c[2]);
      }
    });

    test('direction of nearly antipodal equatorial solution', () {
      final geod = Geodesic.WGS84;
      // [lat1, lat2, azi1, azi2]
      final cases = [
        [0.0, 0.0, 56.0, 124.0],
        [-0.0, -0.0, 124.0, 56.0],
      ];
      for (final c in cases) {
        final inv = geod.inverse(c[0], 0.0, c[1], 179.5);
        expect(inv.azi1, closeTo(c[2], 0.5));
        expect(inv.azi2, closeTo(c[3], 0.5));
      }
    });

    test('direction of the exact antipodal equatorial path', () {
      final geod = Geodesic.WGS84;
      // [lat1, lat2, lon2, expectedAzi1, expectedAzi2]
      final cases = [
        [0.0, 0.0, 180.0, 0.0, 180.0],
        [-0.0, -0.0, 180.0, 180.0, 0.0],
        [0.0, 0.0, -180.0, -0.0, -180.0],
        [-0.0, -0.0, -180.0, -180.0, -0.0],
      ];
      for (final c in cases) {
        final inv = geod.inverse(c[0], 0.0, c[1], c[2]);
        expectExact(inv.azi1, c[3]);
        expectExact(inv.azi2, c[4]);
      }
    });

    test('antipodal points on the equator with prolate ellipsoid', () {
      final geod = Geodesic(6.4e6, -1 / 300);
      // [lon2, azi1/2]
      final cases = [
        [180.0, 90.0],
        [-180.0, -90.0],
      ];
      for (final c in cases) {
        final inv = geod.inverse(0.0, 0.0, 0.0, c[0]);
        expectExact(inv.azi1, c[1]);
        expectExact(inv.azi2, c[1]);
      }
    });

    test('meridional azimuths for the direct problem', () {
      final geod = Geodesic.WGS84;
      // [azi1, expectedLon2, expectedAzi2]
      final cases = [
        [0.0, 180.0, 180.0],
        [-0.0, -180.0, -180.0],
        [180.0, 180.0, 0.0],
        [-180.0, -180.0, -0.0],
      ];
      for (final c in cases) {
        final dir = geod.direct(0.0, 0.0, c[0], 15e6, gLongUnroll);
        expectExact(dir.lon2, c[1]);
        expectExact(dir.azi2, c[2]);
      }
    });
  });
}
