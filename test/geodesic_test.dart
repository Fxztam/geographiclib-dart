// geodesic_test.dart
// Dart port of geodesictest.js from geographiclib-geodesic v2.2.0
// Tests Geodesic, GeodesicLine, and PolygonArea.
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

import 'package:test/test.dart';
import 'package:geographiclib_dart/geodesic.dart';

// ---------------------------------------------------------------------------
// 20 reference test cases: [lat1,lon1,azi1, lat2,lon2,azi2, s12,a12,m12,M12,M21,S12]
// ---------------------------------------------------------------------------
const List<List<double>> _testcases = [
  [35.60777, -139.44815, 111.098748429560326,
   -11.17491, -69.95921, 129.289270889708762,
   8935244.5604818305, 80.50729714281974, 6273170.2055303837,
   0.16606318447386067, 0.16479116945612937, 12841384694976.432],
  [55.52454, 106.05087, 22.020059880982801,
   77.03196, 197.18234, 109.112041110671519,
   4105086.1713924406, 36.892740690445894, 3828869.3344387607,
   0.80076349608092607, 0.80101006984201008, 61674961290615.615],
  [-21.97856, 142.59065, -32.44456876433189,
   41.84138, 98.56635, -41.84359951440466,
   8394328.894657671, 75.62930491011522, 6161154.5773110616,
   0.24816339233950381, 0.24930251203627892, -6637997720646.717],
  [-66.99028, 112.2363, 173.73491240878403,
   -12.70631, 285.90344, 2.512956620913668,
   11150344.2312080241, 100.278634181155759, 6289939.5670446687,
   -0.17199490274700385, -0.17722569526345708, -121287239862139.744],
  [-17.42761, 173.34268, -159.033557661192928,
   -15.84784, 5.93557, -20.787484651536988,
   16076603.1631180673, 144.640108810286253, 3732902.1583877189,
   -0.81273638700070476, -0.81299800519154474, 97825992354058.708],
  [32.84994, 48.28919, 150.492927788121982,
   -56.28556, 202.29132, 48.113449399816759,
   16727068.9438164461, 150.565799985466607, 3147838.1910180939,
   -0.87334918086923126, -0.86505036767110637, -72445258525585.010],
  [6.96833, 52.74123, 92.581585386317712,
   -7.39675, 206.17291, 90.721692165923907,
   17102477.2496958388, 154.147366239113561, 2772035.6169917581,
   -0.89991282520302447, -0.89986892177110739, -1311796973197.995],
  [-50.56724, -16.30485, -105.439679907590164,
   -33.56571, -94.97412, -47.348547835650331,
   6455670.5118668696, 58.083719495371259, 5409150.7979815838,
   0.53053508035997263, 0.52988722644436602, 41071447902810.047],
  [-58.93002, -8.90775, 140.965397902500679,
   -8.91104, 133.13503, 19.255429433416599,
   11756066.0219864627, 105.755691241406877, 6151101.2270708536,
   -0.26548622269867183, -0.27068483874510741, -86143460552774.735],
  [-68.82867, -74.28391, 93.774347763114881,
   -50.63005, -8.36685, 34.65564085411343,
   3956936.926063544, 35.572254987389284, 3708890.9544062657,
   0.81443963736383502, 0.81420859815358342, -41845309450093.787],
  [-10.62672, -32.0898, -86.426713286747751,
   5.883, -134.31681, -80.473780971034875,
   11470869.3864563009, 103.387395634504061, 6184411.6622659713,
   -0.23138683500430237, -0.23155097622286792, 4198803992123.548],
  [-21.76221, 166.90563, 29.319421206936428,
   48.72884, 213.97627, 43.508671946410168,
   9098627.3986554915, 81.963476716121964, 6299240.9166992283,
   0.13965943368590333, 0.14152969707656796, 10024709850277.476],
  [-19.79938, -174.47484, 71.167275780171533,
   -11.99349, -154.35109, 65.589099775199228,
   2319004.8601169389, 20.896611684802389, 2267960.8703918325,
   0.93427001867125849, 0.93424887135032789, -3935477535005.785],
  [-11.95887, -116.94513, 92.712619830452549,
   4.57352, 7.16501, 78.64960934409585,
   13834722.5801401374, 124.688684161089762, 5228093.177931598,
   -0.56879356755666463, -0.56918731952397221, -9919582785894.853],
  [-87.85331, 85.66836, -65.120313040242748,
   66.48646, 16.09921, -4.888658719272296,
   17286615.3147144645, 155.58592449699137, 2635887.4729110181,
   -0.90697975771398578, -0.91095608883042767, 42667211366919.534],
  [1.74708, 128.32011, -101.584843631173858,
   -11.16617, 11.87109, -86.325793296437476,
   12942901.1241347408, 116.650512484301857, 5682744.8413270572,
   -0.44857868222697644, -0.44824490340007729, 10763055294345.653],
  [-25.72959, -144.90758, -153.647468693117198,
   -57.70581, -269.17879, -48.343983158876487,
   9413446.7452453107, 84.664533838404295, 6356176.6898881281,
   0.09492245755254703, 0.09737058264766572, 74515122850712.444],
  [-41.22777, 122.32875, 14.285113402275739,
   -7.57291, 130.37946, 10.805303085187369,
   3812686.035106021, 34.34330804743883, 3588703.8812128856,
   0.82605222593217889, 0.82572158200920196, -2456961531057.857],
  [11.01307, 138.25278, 79.43682622782374,
   6.62726, 247.05981, 103.708090215522657,
   11911190.819018408, 107.341669954114577, 6070904.722786735,
   -0.29767608923657404, -0.29785143390252321, 17121631423099.696],
  [-29.47124, 95.14681, -163.779130441688382,
   -27.46601, -69.15955, -15.909335945554969,
   13487015.8381145492, 121.294026715742277, 5481428.9945736388,
   -0.51527225545373252, -0.51556587964721788, 104679964020340.318],
];

// ---------------------------------------------------------------------------
// Approx helper
// ---------------------------------------------------------------------------
void approx(double actual, double expected, double tol) {
  expect((actual - expected).abs(), lessThanOrEqualTo(tol),
      reason: '$actual ≈ $expected ± $tol');
}

void main() {
  final geod = Geodesic.WGS84;

  // ── Reference table checks ──────────────────────────────────────────────
  group('GeodesicTest', () {
    void checkInverse(List<double> l) {
      final lat1 = l[0], lon1 = l[1], azi1 = l[2];
      final lat2 = l[3], lon2 = l[4], azi2 = l[5];
      final s12 = l[6], a12 = l[7], m12 = l[8];
      final M12 = l[9], M21 = l[10], S12 = l[11];
      final inv = geod.inverse(
          lat1, lon1, lat2, GeoMath.angNormalize(lon2), gAll | gLongUnroll);
      approx(inv.lon2, lon2, 1e-13);
      approx(inv.azi1, azi1, 1e-13);
      approx(inv.azi2, azi2, 1e-13);
      approx(inv.s12, s12, 1e-8);
      approx(inv.a12, a12, 1e-13);
      approx(inv.m12, m12, 1e-8);
      approx(inv.M12, M12, 1e-15);
      approx(inv.M21, M21, 1e-15);
      approx(inv.S12, S12, 0.1);
    }

    void checkDirect(List<double> l) {
      final lat1 = l[0], lon1 = l[1], azi1 = l[2];
      final lat2 = l[3], lon2 = l[4], azi2 = l[5];
      final s12 = l[6], a12 = l[7], m12 = l[8];
      final M12 = l[9], M21 = l[10], S12 = l[11];
      final dir = geod.direct(lat1, lon1, azi1, s12, gAll | gLongUnroll);
      approx(dir.lat2, lat2, 1e-13);
      approx(dir.lon2, lon2, 1e-13);
      approx(dir.azi2, azi2, 1e-13);
      approx(dir.a12, a12, 1e-13);
      approx(dir.m12, m12, 1e-8);
      approx(dir.M12, M12, 1e-15);
      approx(dir.M21, M21, 1e-15);
      approx(dir.S12, S12, 0.1);
    }

    void checkArcDirect(List<double> l) {
      final lat1 = l[0], lon1 = l[1], azi1 = l[2];
      final lat2 = l[3], lon2 = l[4], azi2 = l[5];
      final s12 = l[6], a12 = l[7], m12 = l[8];
      final M12 = l[9], M21 = l[10], S12 = l[11];
      final dir = geod.arcDirect(lat1, lon1, azi1, a12, gAll | gLongUnroll);
      approx(dir.lat2, lat2, 1e-13);
      approx(dir.lon2, lon2, 1e-13);
      approx(dir.azi2, azi2, 1e-13);
      approx(dir.s12, s12, 1e-8);
      approx(dir.m12, m12, 1e-8);
      approx(dir.M12, M12, 1e-15);
      approx(dir.M21, M21, 1e-15);
      approx(dir.S12, S12, 0.1);
    }

    test('check inverse (20 reference cases)', () {
      for (final tc in _testcases) {
        checkInverse(tc);
      }
    });

    test('check direct (20 reference cases)', () {
      for (final tc in _testcases) {
        checkDirect(tc);
      }
    });

    test('check arcDirect (20 reference cases)', () {
      for (final tc in _testcases) {
        checkArcDirect(tc);
      }
    });
  });

  // ── GeodesicSolve tests ─────────────────────────────────────────────────
  group('GeodesicSolve', () {
    test('GeodSolve0', () {
      final inv = geod.inverse(40.6, -73.8, 49.01666667, 2.55);
      approx(inv.azi1, 53.47022, 0.5e-5);
      approx(inv.azi2, 111.59367, 0.5e-5);
      approx(inv.s12, 5853226, 0.5);
    });

    test('GeodSolve1', () {
      final dir = geod.direct(40.63972222, -73.77888889, 53.5, 5850e3);
      approx(dir.lat2, 49.01467, 0.5e-5);
      approx(dir.lon2, 2.56106, 0.5e-5);
      approx(dir.azi2, 111.62947, 0.5e-5);
    });

    test('GeodSolve2', () {
      // Check fix for antipodal prolate bug found 2010-09-04
      final geod2 = Geodesic(6.4e6, -1 / 150);
      var inv = geod2.inverse(0.07476, 0, -0.07476, 180);
      approx(inv.azi1, 90.00078, 0.5e-5);
      approx(inv.azi2, 90.00078, 0.5e-5);
      approx(inv.s12, 20106193, 0.5);
      inv = geod2.inverse(0.1, 0, -0.1, 180);
      approx(inv.azi1, 90.00105, 0.5e-5);
      approx(inv.azi2, 90.00105, 0.5e-5);
      approx(inv.s12, 20106193, 0.5);
    });

    test('GeodSolve4', () {
      // Check fix for short line bug found 2010-05-21
      final inv = geod.inverse(36.493349428792, 0,
          36.49334942879201, 0.0000008);
      approx(inv.s12, 0.072, 0.5e-3);
    });

    test('GeodSolve5', () {
      // Check fix for point2=pole bug found 2010-05-03
      final dir = geod.direct(0.01777745589997, 30, 0, 10e6);
      approx(dir.lat2, 90, 0.5e-5);
      if (dir.lon2 < 0) {
        approx(dir.lon2, -150, 0.5e-5);
        approx(dir.azi2.abs(), 180, 0.5e-5);
      } else {
        approx(dir.lon2, 30, 0.5e-5);
        approx(dir.azi2, 0, 0.5e-5);
      }
    });

    test('GeodSolve6', () {
      // Check fix for volatile sbet12a bug found 2011-06-25
      var inv = geod.inverse(88.202499451857, 0,
          -88.202499451857, 179.981022032992859592);
      approx(inv.s12, 20003898.214, 0.5e-3);
      inv = geod.inverse(89.262080389218, 0,
          -89.262080389218, 179.992207982775375662);
      approx(inv.s12, 20003925.854, 0.5e-3);
      inv = geod.inverse(89.333123580033, 0,
          -89.333123580032997687, 179.99295812360148422);
      approx(inv.s12, 20003926.881, 0.5e-3);
    });

    test('GeodSolve9', () {
      final inv = geod.inverse(56.320923501171, 0,
          -56.320923501171, 179.664747671772880215);
      approx(inv.s12, 19993558.287, 0.5e-3);
    });

    test('GeodSolve10', () {
      final inv = geod.inverse(52.784459512564, 0,
          -52.784459512563990912, 179.634407464943777557);
      approx(inv.s12, 19991596.095, 0.5e-3);
    });

    test('GeodSolve11', () {
      final inv = geod.inverse(48.522876735459, 0,
          -48.52287673545898293, 179.599720456223079643);
      approx(inv.s12, 19989144.774, 0.5e-3);
    });

    test('GeodSolve12', () {
      // Check fix for inverse geodesics on extreme prolate/oblate ellipsoids
      final geod2 = Geodesic(89.8, -1.83);
      final inv = geod2.inverse(0, 0, -10, 160);
      approx(inv.azi1, 120.27, 1e-2);
      approx(inv.azi2, 105.15, 1e-2);
      approx(inv.s12, 266.7, 1e-1);
    });

    test('GeodSolve14', () {
      // Check fix for inverse ignoring lon12 = nan
      final inv = geod.inverse(0, 0, 1, double.nan);
      expect(inv.azi1.isNaN, isTrue);
      expect(inv.azi2.isNaN, isTrue);
      expect(inv.s12.isNaN, isTrue);
    });

    test('GeodSolve15', () {
      // Check fix for Math::eatanhe for e^2 < 0
      final geod2 = Geodesic(6.4e6, -1 / 150);
      final dir = geod2.direct(1, 2, 3, 4, gArea);
      approx(dir.S12, 23700, 0.5);
    });

    test('GeodSolve17', () {
      // Check fix for LONG_UNROLL bug found on 2015-05-07
      var dir = geod.direct(40, -75, -10, 2e7, gLongUnroll);
      approx(dir.lat2, -39, 1);
      approx(dir.lon2, -254, 1);
      approx(dir.azi2, -170, 1);
      final line = geod.line(40, -75, -10);
      dir = line.position(2e7, gLongUnroll);
      approx(dir.lat2, -39, 1);
      approx(dir.lon2, -254, 1);
      approx(dir.azi2, -170, 1);
      dir = geod.direct(40, -75, -10, 2e7);
      approx(dir.lat2, -39, 1);
      approx(dir.lon2, 105, 1);
      approx(dir.azi2, -170, 1);
      dir = line.position(2e7);
      approx(dir.lat2, -39, 1);
      approx(dir.lon2, 105, 1);
      approx(dir.azi2, -170, 1);
    });

    test('GeodSolve26', () {
      // Check 0/0 problem with area calculation on sphere 2015-09-08
      final geod2 = Geodesic(6.4e6, 0);
      final inv = geod2.inverse(1, 2, 3, 4, gArea);
      approx(inv.S12, 49911046115, 0.5);
    });

    test('GeodSolve28', () {
      final geod2 = Geodesic(6.4e6, 0.1);
      final dir = geod2.direct(1, 2, 10, 5e6);
      approx(dir.a12, 48.55570690, 0.5e-8);
    });

    test('GeodSolve29', () {
      var dir = geod.inverse(0, 539, 0, 181);
      approx(dir.lon1, 179, 1e-10);
      approx(dir.lon2, -179, 1e-10);
      approx(dir.s12, 222639, 0.5);
      dir = geod.inverse(0, 539, 0, 181, gLongUnroll);
      approx(dir.lon1, 539, 1e-10);
      approx(dir.lon2, 541, 1e-10);
      approx(dir.s12, 222639, 0.5);
    });

    test('GeodSolve33', () {
      var inv = geod.inverse(0, 0, 0, 179);
      approx(inv.azi1, 90.00000, 0.5e-5);
      approx(inv.azi2, 90.00000, 0.5e-5);
      approx(inv.s12, 19926189, 0.5);
      inv = geod.inverse(0, 0, 0, 179.5);
      approx(inv.azi1, 55.96650, 0.5e-5);
      approx(inv.azi2, 124.03350, 0.5e-5);
      approx(inv.s12, 19980862, 0.5);
      inv = geod.inverse(0, 0, 0, 180);
      approx(inv.azi1, 0.00000, 0.5e-5);
      approx(inv.azi2.abs(), 180.00000, 0.5e-5);
      approx(inv.s12, 20003931, 0.5);
      inv = geod.inverse(0, 0, 1, 180);
      approx(inv.azi1, 0.00000, 0.5e-5);
      approx(inv.azi2.abs(), 180.00000, 0.5e-5);
      approx(inv.s12, 19893357, 0.5);

      final sph = Geodesic(6.4e6, 0);
      inv = sph.inverse(0, 0, 0, 179);
      approx(inv.azi1, 90.00000, 0.5e-5);
      approx(inv.azi2, 90.00000, 0.5e-5);
      approx(inv.s12, 19994492, 0.5);
      inv = sph.inverse(0, 0, 0, 180);
      approx(inv.azi1, 0.00000, 0.5e-5);
      approx(inv.azi2.abs(), 180.00000, 0.5e-5);
      approx(inv.s12, 20106193, 0.5);
      inv = sph.inverse(0, 0, 1, 180);
      approx(inv.azi1, 0.00000, 0.5e-5);
      approx(inv.azi2.abs(), 180.00000, 0.5e-5);
      approx(inv.s12, 19994492, 0.5);

      final pro = Geodesic(6.4e6, -1 / 300);
      inv = pro.inverse(0, 0, 0, 179);
      approx(inv.azi1, 90.00000, 0.5e-5);
      approx(inv.azi2, 90.00000, 0.5e-5);
      approx(inv.s12, 19994492, 0.5);
      inv = pro.inverse(0, 0, 0, 180);
      approx(inv.azi1, 90.00000, 0.5e-5);
      approx(inv.azi2, 90.00000, 0.5e-5);
      approx(inv.s12, 20106193, 0.5);
      inv = pro.inverse(0, 0, 0.5, 180);
      approx(inv.azi1, 33.02493, 0.5e-5);
      approx(inv.azi2, 146.97364, 0.5e-5);
      approx(inv.s12, 20082617, 0.5);
      inv = pro.inverse(0, 0, 1, 180);
      approx(inv.azi1, 0.00000, 0.5e-5);
      approx(inv.azi2.abs(), 180.00000, 0.5e-5);
      approx(inv.s12, 20027270, 0.5);
    });

    test('GeodSolve55', () {
      var inv = geod.inverse(double.nan, 0, 0, 90);
      expect(inv.azi1.isNaN, isTrue);
      expect(inv.azi2.isNaN, isTrue);
      expect(inv.s12.isNaN, isTrue);
      inv = geod.inverse(double.nan, 0, 90, 9);
      expect(inv.azi1.isNaN, isTrue);
      expect(inv.azi2.isNaN, isTrue);
      expect(inv.s12.isNaN, isTrue);
    });

    test('GeodSolve59', () {
      final inv = geod.inverse(5, 0.00000000000001, 10, 180);
      approx(inv.azi1, 0.000000000000035, 1.5e-14);
      approx(inv.azi2, 179.99999999999996, 1.5e-14);
      approx(inv.s12, 18345191.174332713, 5e-9);
    });

    test('GeodSolve61', () {
      // Make sure small negative azimuths are west-going
      var dir = geod.direct(45, 0, -0.000000000000000003, 1e7, gLongUnroll);
      approx(dir.lat2, 45.30632, 0.5e-5);
      approx(dir.lon2, -180, 0.5e-5);
      approx(dir.azi2.abs(), 180, 0.5e-5);

      final line = geod.inverseLine(45, 0, 80, -0.000000000000000003);
      dir = line.position(1e7, gLongUnroll);
      approx(dir.lat2, 45.30632, 0.5e-5);
      approx(dir.lon2, -180, 0.5e-5);
      approx(dir.azi2.abs(), 180, 0.5e-5);
    });

    test('GeodSolve65', () {
      final line = geod.inverseLine(
          30, -0.000000000000000001, -31, 180, gAll);
      var dir = line.position(1e7, gAll | gLongUnroll);
      approx(dir.lat1, 30.00000, 0.5e-5);
      approx(dir.lon1, -0.00000, 0.5e-5);
      approx(dir.azi1.abs(), 180.00000, 0.5e-5);
      approx(dir.lat2, -60.23169, 0.5e-5);
      approx(dir.lon2, -0.00000, 0.5e-5);
      approx(dir.azi2.abs(), 180.00000, 0.5e-5);
      approx(dir.s12, 10000000, 0.5);
      approx(dir.a12, 90.06544, 0.5e-5);
      approx(dir.m12, 6363636, 0.5);
      approx(dir.M12, -0.0012834, 0.5e-7);
      approx(dir.M21, 0.0013749, 0.5e-7);
      approx(dir.S12, 0, 0.5);

      dir = line.position(2e7, gAll | gLongUnroll);
      approx(dir.lat1, 30.00000, 0.5e-5);
      approx(dir.lon1, -0.00000, 0.5e-5);
      approx(dir.azi1.abs(), 180.00000, 0.5e-5);
      approx(dir.lat2, -30.03547, 0.5e-5);
      approx(dir.lon2, -180.00000, 0.5e-5);
      approx(dir.azi2, -0.00000, 0.5e-5);
      approx(dir.s12, 20000000, 0.5);
      approx(dir.a12, 179.96459, 0.5e-5);
      approx(dir.m12, 54342, 0.5);
      approx(dir.M12, -1.0045592, 0.5e-7);
      approx(dir.M21, -0.9954339, 0.5e-7);
      approx(dir.S12, 127516405431022, 0.5);
    });

    test('GeodSolve69', () {
      final line = geod.inverseLine(-5, -0.000000000000002, -10, 180);
      var dir = line.position(2e7, gLongUnroll);
      approx(dir.lat2, 4.96445, 0.5e-5);
      approx(dir.lon2, -180.00000, 0.5e-5);
      approx(dir.azi2, -0.00000, 0.5e-5);
      dir = line.position(0.5 * line.s13, gLongUnroll);
      approx(dir.lat2, -87.52461, 0.5e-5);
      approx(dir.lon2, -0.00000, 0.5e-5);
      approx(dir.azi2, -180.00000, 0.5e-5);
    });

    test('GeodSolve71', () {
      // Check that DirectLine sets s13.
      final line = geod.directLine(1, 2, 45, 1e7);
      final dir = line.position(0.5 * line.s13, gLongUnroll);
      approx(dir.lat2, 30.92625, 0.5e-5);
      approx(dir.lon2, 37.54640, 0.5e-5);
      approx(dir.azi2, 55.43104, 0.5e-5);
    });

    test('GeodSolve73', () {
      // Check for backwards from the pole bug reported 2016-02-13
      final dir = geod.direct(90, 10, 180, -1e6);
      approx(dir.lat2, 81.04623, 0.5e-5);
      approx(dir.lon2, -170, 0.5e-5);
      approx(dir.azi2, 0, 0.5e-5);
      expect(GeoMath.copysign(1, dir.azi2), greaterThan(0));
    });

    test('GeodSolve74', () {
      // Check fix for inaccurate areas, bug introduced in v1.46
      final inv = geod.inverse(54.1589, 15.3872, 54.1591, 15.3877, gAll);
      approx(inv.azi1, 55.723110355, 5e-9);
      approx(inv.azi2, 55.723515675, 5e-9);
      approx(inv.s12, 39.527686385, 5e-9);
      approx(inv.a12, 0.000355495, 5e-9);
      approx(inv.m12, 39.527686385, 5e-9);
      approx(inv.M12, 0.999999995, 5e-9);
      approx(inv.M21, 0.999999995, 5e-9);
      approx(inv.S12, 286698586.30197, 5e-4);
    });

    test('GeodSolve76', () {
      // The distance from Wellington and Salamanca (classic Vincenty failure)
      final inv = geod.inverse(
          -(41 + 19 / 60), 174 + 49 / 60, 40 + 58 / 60, -(5 + 30 / 60));
      approx(inv.azi1, 160.39137649664, 0.5e-11);
      approx(inv.azi2, 19.50042925176, 0.5e-11);
      approx(inv.s12, 19960543.857179, 0.5e-6);
    });

    test('GeodSolve78', () {
      // An example where the NGS calculator fails to converge
      final inv = geod.inverse(27.2, 0, -27.1, 179.5);
      approx(inv.azi1, 45.82468716758, 0.5e-11);
      approx(inv.azi2, 134.22776532670, 0.5e-11);
      approx(inv.s12, 19974354.765767, 0.5e-6);
    });

    test('GeodSolve80', () {
      // Computing scale in special cases + zero length geodesic
      var inv = geod.inverse(0, 0, 0, 90, gGeodesicScale);
      approx(inv.M12, -0.00528427534, 0.5e-10);
      approx(inv.M21, -0.00528427534, 0.5e-10);

      inv = geod.inverse(0, 0, 1e-6, 1e-6, gGeodesicScale);
      approx(inv.M12, 1, 0.5e-10);
      approx(inv.M21, 1, 0.5e-10);

      inv = geod.inverse(20.001, 0, 20.001, 0, gAll);
      approx(inv.a12, 0, 1e-13);
      approx(inv.s12, 0, 1e-8);
      approx(inv.azi1, 180, 1e-13);
      approx(inv.azi2, 180, 1e-13);
      approx(inv.m12, 0, 1e-8);
      approx(inv.M12, 1, 1e-15);
      approx(inv.M21, 1, 1e-15);
      approx(inv.S12, 0, 1e-10);
      expect(GeoMath.copysign(1, inv.a12), greaterThan(0));
      expect(GeoMath.copysign(1, inv.s12), greaterThan(0));
      expect(GeoMath.copysign(1, inv.m12), greaterThan(0));

      inv = geod.inverse(90, 0, 90, 180, gAll);
      approx(inv.a12, 0, 1e-13);
      approx(inv.s12, 0, 1e-8);
      approx(inv.azi1, 0, 1e-13);
      approx(inv.azi2, 180, 1e-13);
      approx(inv.m12, 0, 1e-8);
      approx(inv.M12, 1, 1e-15);
      approx(inv.M21, 1, 1e-15);
      approx(inv.S12, 127516405431022, 0.5);

      // Incapable line that can't take distance as input
      final ln = geod.line(1, 2, 90, gLatitude);
      final dir = ln.position(1000, gNone);
      expect(dir.a12.isNaN, isTrue);
    });

    test('GeodSolve84', () {
      // Tests for range errors with {fmod,sin,cos}(inf)
      var dir = geod.direct(0, 0, 90, double.infinity);
      expect(dir.lat2.isNaN, isTrue);
      expect(dir.lon2.isNaN, isTrue);
      expect(dir.azi2.isNaN, isTrue);
      dir = geod.direct(0, 0, 90, double.nan);
      expect(dir.lat2.isNaN, isTrue);
      expect(dir.lon2.isNaN, isTrue);
      expect(dir.azi2.isNaN, isTrue);
      dir = geod.direct(0, 0, double.infinity, 1000);
      expect(dir.lat2.isNaN, isTrue);
      expect(dir.lon2.isNaN, isTrue);
      expect(dir.azi2.isNaN, isTrue);
      dir = geod.direct(0, 0, double.nan, 1000);
      expect(dir.lat2.isNaN, isTrue);
      expect(dir.lon2.isNaN, isTrue);
      expect(dir.azi2.isNaN, isTrue);
      dir = geod.direct(0, double.infinity, 90, 1000);
      expect(dir.lat2, equals(0));
      expect(dir.lon2.isNaN, isTrue);
      expect(dir.azi2, equals(90));
      dir = geod.direct(0, double.nan, 90, 1000);
      expect(dir.lat2, equals(0));
      expect(dir.lon2.isNaN, isTrue);
      expect(dir.azi2, equals(90));
      dir = geod.direct(double.infinity, 0, 90, 1000);
      expect(dir.lat2.isNaN, isTrue);
      expect(dir.lon2.isNaN, isTrue);
      expect(dir.azi2.isNaN, isTrue);
      dir = geod.direct(double.nan, 0, 90, 1000);
      expect(dir.lat2.isNaN, isTrue);
      expect(dir.lon2.isNaN, isTrue);
      expect(dir.azi2.isNaN, isTrue);
    });

    test('GeodSolve92', () {
      // Check fix for inaccurate hypot with python 3.[89]
      final inv = geod.inverse(
          37.757540000000006, -122.47018, 37.75754, -122.470177);
      approx(inv.azi1, 89.99999923, 1e-7);
      approx(inv.azi2, 90.00000106, 1e-7);
      approx(inv.s12, 0.264, 0.5e-3);
    });

    test('GeodSolve94', () {
      // Check fix for lat2 = nan being treated as lat2 = 0
      final inv = geod.inverse(0, 0, double.nan, 90);
      expect(inv.azi1.isNaN, isTrue);
      expect(inv.azi2.isNaN, isTrue);
      expect(inv.s12.isNaN, isTrue);
    });

    test('GeodSolve96', () {
      // Failure with long doubles (Nowak + Nowak Da Costa 2022)
      final geod2 = Geodesic(6378137, 1 / 298.257222101);
      final inv = geod2.inverse(
          0, 0, 60.0832522871723, 89.8492185074635, gArea);
      approx(inv.S12, 42426932221845, 0.5);
    });

    test('GeodSolve99', () {
      // sincosd(+/-45) inconsistency due to directed rounding
      final inv = geod.inverse(45, 0, -45, 179.572719);
      approx(inv.azi1, 90.00000028, 1e-8);
      approx(inv.azi2, 90.00000028, 1e-8);
      approx(inv.s12, 19987083.007, 0.5e-3);
    });

    test('GeodSolve100', () {
      // Check fix for meridional failure for a strongly prolate ellipsoid
      final geod2 = Geodesic(1e6, -3);
      final inv = geod2.inverse(30, 0, 30, 180);
      approx(inv.azi1, 22.368806, 1.0);
      approx(inv.azi2, 157.631194, 1.0);
      approx(inv.s12, 1074081.6, 1e3);
    });
  });

  // ── Planimeter tests ────────────────────────────────────────────────────
  group('Planimeter', () {
    late PolygonArea polygon;
    late PolygonArea polyline;

    setUp(() {
      polygon = geod.polygon(false);
      polyline = geod.polygon(true);
    });

    PolygonResult planimeter(List<List<double>> points) {
      polygon.clear();
      for (final p in points) {
        polygon.addPoint(p[0], p[1]);
      }
      return polygon.compute(reverse: false, sign: true);
    }

    PolygonResult polyLength(List<List<double>> points) {
      polyline.clear();
      for (final p in points) {
        polyline.addPoint(p[0], p[1]);
      }
      return polyline.compute(reverse: false, sign: true);
    }

    test('Planimeter0', () {
      // Check fix for pole-encircling bug found 2011-03-16
      var a = planimeter([[89, 0], [89, 90], [89, 180], [89, 270]]);
      approx(a.perimeter, 631819.8745, 1e-4);
      approx(a.area, 24952305678, 1);

      a = planimeter([[-89, 0], [-89, 90], [-89, 180], [-89, 270]]);
      approx(a.perimeter, 631819.8745, 1e-4);
      approx(a.area, -24952305678, 1);

      a = planimeter([[0, -1], [-1, 0], [0, 1], [1, 0]]);
      approx(a.perimeter, 627598.2731, 1e-4);
      approx(a.area, 24619419146, 1);

      a = planimeter([[90, 0], [0, 0], [0, 90]]);
      approx(a.perimeter, 30022685, 1);
      approx(a.area, 63758202715511, 1);

      final pl = polyLength([[90, 0], [0, 0], [0, 90]]);
      approx(pl.perimeter, 20020719, 1);
      expect(pl.area.isNaN, isTrue);
    });

    test('Planimeter5', () {
      // Check fix for Planimeter pole crossing bug found 2011-06-24
      final a = planimeter([[89, 0.1], [89, 90.1], [89, -179.9]]);
      approx(a.perimeter, 539297, 1);
      approx(a.area, 12476152838.5, 1);
    });

    test('Planimeter6', () {
      // Check fix for Planimeter lon12 rounding bug found 2012-12-03
      var a = planimeter([[9, -0.00000000000001], [9, 180], [9, 0]]);
      approx(a.perimeter, 36026861, 1);
      approx(a.area, 0, 1);
      a = planimeter([[9, 0.00000000000001], [9, 0], [9, 180]]);
      approx(a.perimeter, 36026861, 1);
      approx(a.area, 0, 1);
      a = planimeter([[9, 0.00000000000001], [9, 180], [9, 0]]);
      approx(a.perimeter, 36026861, 1);
      approx(a.area, 0, 1);
      a = planimeter([[9, -0.00000000000001], [9, 0], [9, 180]]);
      approx(a.perimeter, 36026861, 1);
      approx(a.area, 0, 1);
    });

    test('Planimeter12', () {
      // Area of arctic circle
      final a = planimeter(
          [[66.562222222, 0], [66.562222222, 180], [66.562222222, 360]]);
      approx(a.perimeter, 10465729, 1);
      approx(a.area, 0, 1);
    });

    test('Planimeter12r', () {
      // Reverse area of arctic circle
      final a = planimeter(
          [[66.562222222, -0.0], [66.562222222, -180], [66.562222222, -360]]);
      approx(a.perimeter, 10465729, 1);
      approx(a.area, 0, 1);
    });

    test('Planimeter13', () {
      // Check encircling pole twice
      final a = planimeter([
        [89, -360], [89, -240], [89, -120],
        [89, 0], [89, 120], [89, 240]
      ]);
      approx(a.perimeter, 1160741, 1);
      approx(a.area, 32415230256, 1);
    });

    test('Planimeter15', () {
      // Coverage tests for reverse/sign combinations
      final lat = [2.0, 1.0, 3.0];
      final lon = [1.0, 2.0, 3.0];
      const r = 18454562325.45119;
      const a0 = 510065621724088.5093;

      polygon.clear();
      polygon.addPoint(lat[0], lon[0]);
      polygon.addPoint(lat[1], lon[1]);

      var a = polygon.testPoint(lat[2], lon[2], reverse: false, sign: true);
      approx(a.area, r, 0.5);
      a = polygon.testPoint(lat[2], lon[2], reverse: false, sign: false);
      approx(a.area, r, 0.5);
      a = polygon.testPoint(lat[2], lon[2], reverse: true, sign: true);
      approx(a.area, -r, 0.5);
      a = polygon.testPoint(lat[2], lon[2], reverse: true, sign: false);
      approx(a.area, a0 - r, 0.5);

      final inv = geod.inverse(lat[1], lon[1], lat[2], lon[2]);
      a = polygon.testEdge(inv.azi1, inv.s12, reverse: false, sign: true);
      approx(a.area, r, 0.5);
      a = polygon.testEdge(inv.azi1, inv.s12, reverse: false, sign: false);
      approx(a.area, r, 0.5);
      a = polygon.testEdge(inv.azi1, inv.s12, reverse: true, sign: true);
      approx(a.area, -r, 0.5);
      a = polygon.testEdge(inv.azi1, inv.s12, reverse: true, sign: false);
      approx(a.area, a0 - r, 0.5);

      polygon.addPoint(lat[2], lon[2]);
      a = polygon.compute(reverse: false, sign: true);
      approx(a.area, r, 0.5);
      a = polygon.compute(reverse: false, sign: false);
      approx(a.area, r, 0.5);
      a = polygon.compute(reverse: true, sign: true);
      approx(a.area, -r, 0.5);
      a = polygon.compute(reverse: true, sign: false);
      approx(a.area, a0 - r, 0.5);
    });

    test('Planimeter19', () {
      // Degenerate polygon coverage tests
      polygon.clear();
      var a = polygon.compute(reverse: false, sign: true);
      expect(a.area, equals(0));
      expect(a.perimeter, equals(0));

      a = polygon.testPoint(1, 1, reverse: false, sign: true);
      expect(a.area, equals(0));
      expect(a.perimeter, equals(0));

      a = polygon.testEdge(90, 1000, reverse: false, sign: true);
      expect(a.area.isNaN, isTrue);
      expect(a.perimeter.isNaN, isTrue);

      polygon.addPoint(1, 1);
      a = polygon.compute(reverse: false, sign: true);
      expect(a.area, equals(0));
      expect(a.perimeter, equals(0));

      polyline.clear();
      var pl = polyline.compute(reverse: false, sign: true);
      expect(pl.perimeter, equals(0));

      pl = polyline.testPoint(1, 1, reverse: false, sign: true);
      expect(pl.perimeter, equals(0));

      pl = polyline.testEdge(90, 1000, reverse: false, sign: true);
      expect(pl.perimeter.isNaN, isTrue);

      polyline.addPoint(1, 1);
      pl = polyline.compute(reverse: false, sign: true);
      expect(pl.perimeter, equals(0));

      polygon.addPoint(1, 1);
      pl = polyline.testEdge(90, 1000, reverse: false, sign: true);
      approx(pl.perimeter, 1000, 1e-10);

      pl = polyline.testPoint(2, 2, reverse: false, sign: true);
      approx(pl.perimeter, 156876.149, 0.5e-3);
    });

    test('Planimeter21', () {
      // Multiple circlings of pole
      const lat = 45.0;
      const azi = 39.2144607176828184218;
      const s = 8420705.40957178156285;
      const r = 39433884866571.4277;
      const a0 = 510065621724088.5093;

      polygon.clear();
      polygon.addPoint(lat, 60);
      polygon.addPoint(lat, 180);
      polygon.addPoint(lat, -60);
      polygon.addPoint(lat, 60);
      polygon.addPoint(lat, 180);
      polygon.addPoint(lat, -60);

      for (int i = 3; i <= 4; ++i) {
        polygon.addPoint(lat, 60);
        polygon.addPoint(lat, 180);

        var a = polygon.testPoint(lat, -60, reverse: false, sign: true);
        approx(a.area, i * r, 0.5);
        a = polygon.testPoint(lat, -60, reverse: false, sign: false);
        approx(a.area, i * r, 0.5);
        a = polygon.testPoint(lat, -60, reverse: true, sign: true);
        approx(a.area, -i * r, 0.5);
        a = polygon.testPoint(lat, -60, reverse: true, sign: false);
        approx(a.area, -i * r + a0, 0.5);

        a = polygon.testEdge(azi, s, reverse: false, sign: true);
        approx(a.area, i * r, 0.5);
        a = polygon.testEdge(azi, s, reverse: false, sign: false);
        approx(a.area, i * r, 0.5);
        a = polygon.testEdge(azi, s, reverse: true, sign: true);
        approx(a.area, -i * r, 0.5);
        a = polygon.testEdge(azi, s, reverse: true, sign: false);
        approx(a.area, -i * r + a0, 0.5);

        polygon.addPoint(lat, -60);
        a = polygon.compute(reverse: false, sign: true);
        approx(a.area, i * r, 0.5);
        a = polygon.compute(reverse: false, sign: false);
        approx(a.area, i * r, 0.5);
        a = polygon.compute(reverse: true, sign: true);
        approx(a.area, -i * r, 0.5);
        a = polygon.compute(reverse: true, sign: false);
        approx(a.area, -i * r + a0, 0.5);
      }
    });

    test('Planimeter29', () {
      // Check fix to transitdirect vs transit zero handling inconsistency
      polygon.clear();
      polygon.addPoint(0, 0);
      polygon.addEdge(90, 1000);
      polygon.addEdge(0, 1000);
      polygon.addEdge(-90, 1000);
      final a = polygon.compute(reverse: false, sign: true);
      approx(a.area, 1000000, 0.01);
    });

    test('check TestEdge', () {
      // Check fix of bugs found by threepointone, 2015-10-14
      polygon.clear();
      polygon.addPoint(33, 44);
      polygon.testEdge(90, 10e3, reverse: false, sign: true);
      polygon.addEdge(90, 10e3);
    });
  });
}
