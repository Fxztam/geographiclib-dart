// geodtest_cpp.dart
// Dart port of C++ test cases from:
//   https://github.com/geographiclib/geographiclib/blob/master/tests/CMakeLists.txt
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
// Scope: GeodSolve and Planimeter tests that are defined in the C++ test
// suite but are NOT present in geodesic_test.dart (the port of geodesictest.js).
//
// Why tests are absent from the JS port:
//   a) The JS combines several C++ sub-tests into one numbered test.
//      (e.g. C++ GeodSolve6/7/8 are all in JS GeodSolve6.)
//   b) GeodesicExact variants (C++ -E flag) — Dart has no GeodesicExact.
//   c) Rhumb-line Planimeter tests (C++ -R flag) — Dart has no Rhumb.
//   d) Command-line DMS-parsing tests (C++ GeodSolve25) — not an API test.
//
// Every test here is labelled with its C++ test name for traceability.

import 'package:test/test.dart';
import 'package:geographiclib_dart/geodesic.dart';

void main() {
  final geod = Geodesic.WGS84;

  // ── GeodSolve tests ──────────────────────────────────────────────────────
  group('GeodSolve (C++ additions)', () {
    // GeodSolve3 — Another check for antipodal prolate bug found 2010-09-04.
    // The JS GeodSolve2 already tests input "0.07476 0 -0.07476 180" and
    // "0.1 0 -0.1 180" in one combined test; C++ splits them as GeodSolve2
    // and GeodSolve3.
    test('GeodSolve3', () {
      final geod2 = Geodesic(6.4e6, -1 / 150.0);
      final inv = geod2.inverse(0.1, 0, -0.1, 180);
      expect(inv.azi1, closeTo(90.00105, 0.5e-5));
      expect(inv.azi2, closeTo(90.00105, 0.5e-5));
      expect(inv.s12,  closeTo(20106193, 0.5));
    });

    // GeodSolve7 — Check fix for volatile sbet12a bug found 2011-06-25
    // (gcc 4.4.4 x86 -O3).  Found again on 2012-03-27 with tdm-mingw32
    // (g++ 4.6.1).  Second sub-case; first is in JS GeodSolve6.
    test('GeodSolve7', () {
      final inv = geod.inverse(
          89.262080389218, 0, -89.262080389218, 179.992207982775375662);
      expect(inv.s12, closeTo(20003925.854, 0.5e-3));
    });

    // GeodSolve8 — Third sub-case of the same volatile sbet12a fix.
    test('GeodSolve8', () {
      final inv = geod.inverse(
          89.333123580033, 0, -89.333123580032997687, 179.99295812360148422);
      expect(inv.s12, closeTo(20003926.881, 0.5e-3));
    });
  });

  // ── Planimeter tests ─────────────────────────────────────────────────────
  group('Planimeter (C++ additions)', () {
    late PolygonArea polygon;
    late PolygonArea polyline;

    setUp(() {
      polygon  = geod.polygon(false);
      polyline = geod.polygon(true);
    });

    // Helper: compute polygon area with default sign conventions.
    PolygonResult planimeter(List<List<double>> pts,
        {bool reverse = false, bool sign = true}) {
      polygon.clear();
      for (final p in pts) polygon.addPoint(p[0], p[1]);
      return polygon.compute(reverse: reverse, sign: sign);
    }

    // Helper: compute polyline length.
    PolygonResult polyLength(List<List<double>> pts) {
      polyline.clear();
      for (final p in pts) polyline.addPoint(p[0], p[1]);
      return polyline.compute(reverse: false, sign: true);
    }

    // Planimeter1 — Check fix for pole-encircling bug (south variant with
    // reversed orientation).
    // C++ test: `Planimeter -r --input-string "-89 0;-89 90;-89 180;-89 270"`
    // The JS Planimeter0 tests this polygon with Compute(false, true) giving
    // area = −24952305678; the C++ Planimeter1 uses `-r` (reverse=true),
    // which flips the sign to +24952305678.
    test('Planimeter1', () {
      final pts = [[-89.0, 0.0], [-89.0, 90.0], [-89.0, 180.0], [-89.0, 270.0]];
      // reverse=true (C++ -r flag) → sign is flipped relative to Planimeter0
      final a = planimeter(pts, reverse: true, sign: true);
      expect(a.perimeter, closeTo(631819.8745, 1e-4));
      expect(a.area,      closeTo(24952305678, 1.0));
    });

    // Planimeter2 — Diamond polygon.
    // C++ test: `Planimeter --input-string "0 -1;-1 0;0 1;1 0"`
    // Sub-case of JS Planimeter0; C++ lists it as a standalone test.
    test('Planimeter2', () {
      final a = planimeter([[0, -1], [-1, 0], [0, 1], [1, 0]]);
      expect(a.perimeter, closeTo(627598.2731, 1e-4));
      expect(a.area,      closeTo(24619419146, 1.0));
    });

    // Planimeter3 — Triangle with north pole vertex.
    // C++ test: `Planimeter --input-string "90 0; 0 0; 0 90"`
    test('Planimeter3', () {
      final a = planimeter([[90, 0], [0, 0], [0, 90]]);
      expect(a.perimeter, closeTo(30022685, 1.0));
      expect(a.area,      closeTo(63758202715511.0, 1.0));
    });

    // Planimeter4 — Same triangle as Planimeter3, but as a polyline.
    // C++ test: `Planimeter -l --input-string "90 0; 0 0; 0 90"`
    test('Planimeter4', () {
      final a = polyLength([[90, 0], [0, 0], [0, 90]]);
      expect(a.perimeter, closeTo(20020719, 1.0));
      expect(a.area.isNaN, isTrue);
    });

    // Planimeter7 — Check fix for Planimeter lon12 rounding bug 2012-12-03.
    // C++ test: second sub-case (positive epsilon, reversed point order).
    // Sub-case of JS Planimeter6; C++ lists all four variations separately.
    test('Planimeter7', () {
      final a = planimeter([[9, 0.00000000000001], [9, 0], [9, 180]]);
      expect(a.perimeter, closeTo(36026861, 1.0));
      expect(a.area,      closeTo(0.0, 1.0));
    });

    // Planimeter8 — Third sub-case.
    test('Planimeter8', () {
      final a = planimeter([[9, 0.00000000000001], [9, 180], [9, 0]]);
      expect(a.perimeter, closeTo(36026861, 1.0));
      expect(a.area,      closeTo(0.0, 1.0));
    });

    // Planimeter9 — Fourth sub-case.
    test('Planimeter9', () {
      final a = planimeter([[9, -0.00000000000001], [9, 0], [9, 180]]);
      expect(a.perimeter, closeTo(36026861, 1.0));
      expect(a.area,      closeTo(0.0, 1.0));
    });

    // ── Tests NOT ported (not applicable to this package) ──
    //
    // Planimeter10 uses C++ flag -R (Rhumb-line polygon) and tests the area
    // of Wyoming; no Rhumb class in this Dart package.
    //
    // Planimeter11 / Planimeter11r also use -R (Rhumb).
    //
    // Planimeter14 uses -w (lon,lat input order swap), identical result to
    // Planimeter13 and adds no algorithmic coverage.
    //
    // GeodSolve13, 16, 18, 27, 30, 34-54, 56, 57, 58, 60, 62-64, 66, 75,
    // 77, 79, 93, 95, 97, 98, 101 all require the C++ GeodesicExact class
    // (invoked via -E flag) which is not implemented in this package.
    //
    // GeodSolve25 tests DMS-format string parsing in the GeodSolve CLI and
    // is not an API-level test.
  });
}
