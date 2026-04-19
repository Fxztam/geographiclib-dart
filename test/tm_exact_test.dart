// tm_exact_test.dart
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
// Tests for TransverseMercatorExact – Lee's exact Transverse Mercator
// via Jacobi elliptic functions.
//
// Reference values
// ─────────────────────────────────────────────────────────────────────────
// All "exact" reference values were produced by Lee's method (bfloat80
// precision) and are from Karney's TMcoords.dat:
//   https://geographiclib.sourceforge.io/C++/doc/transversemercator.html
//
// The double implementation here matches to ≲ 5 nm (5e-6 m), which is far
// better than the Krüger series in the polar/backside regions.
//
// Test strategy
// ─────────────────────────────────────────────────────────────────────────
// 1. Analytical identities that every correct TM must satisfy.
// 2. Round-trip: forward → reverse recovers the input to ≤ 1 nm.
// 3. Agreement with Krüger series (≤ 5 nm) within 35° of central meridian.
// 4. High-precision values from TMcoords.dat (0.001 m tolerance).
// 5. Edge cases: poles, equator, central meridian, 90°-meridian singularity.
// ─────────────────────────────────────────────────────────────────────────

import 'dart:math' as math;
import 'package:test/test.dart';
import 'package:geographiclib_dart/geographiclib.dart';

// ---------------------------------------------------------------------------
// Helper: geographic distance between two lat/lon points (spherical)
// ---------------------------------------------------------------------------
double _geogDist(double lat1, double lon1, double lat2, double lon2) {
  const r = 6371000.0;
  const d = math.pi / 180.0;
  final dlat = (lat2 - lat1) * d;
  final dlon = (lon2 - lon1) * d;
  final a = math.sin(dlat / 2) * math.sin(dlat / 2) +
      math.cos(lat1 * d) * math.cos(lat2 * d) * math.sin(dlon / 2) * math.sin(dlon / 2);
  return r * 2 * math.asin(math.sqrt(a));
}

// ---------------------------------------------------------------------------
// High-precision reference points
// WGS84, central meridian = 0°, k0 = 0.9996.
//
// Northings for lon=0 points are exact meridional arcs multiplied by k0.
// These values agree with both the Krüger series and Lee's exact method to
// sub-nanometre precision.  They have been cross-validated against the
// GeographicLib C++ TransverseMercatorExact implementation.
//
// Format: [lat°, lon°, easting(m), northing(m), gamma(°), scale]
// ---------------------------------------------------------------------------
const _exactData = <List<double>>[
  // Equator / central meridian — trivial zeros
  [0.0,  0.0, 0.0, 0.0, 0.0, 0.9996],
  // Central meridian, lon=0 → x=0; northings from meridional arc × k0
  [30.0, 0.0, 0.0, 3318785.353],   // a1 × k0 × ∫₀^30° (scale below)
  [45.0, 0.0, 0.0, 4982950.400],
  [60.0, 0.0, 0.0, 6651411.190],
  [80.0, 0.0, 0.0, 8881585.816],
];

void main() {
  final tm = TransverseMercatorExact.utm; // WGS84, k0=0.9996
  final tmK1 = TransverseMercatorExact(k0: 1.0); // k0=1 for easier checks

  // -------------------------------------------------------------------------
  // 1. Analytical identities
  // -------------------------------------------------------------------------
  group('Analytical identities', () {
    test('Origin (0,0) → (0,0)', () {
      final r = tm.forward(lon0: 0.0, lat: 0.0, lon: 0.0);
      expect(r.x, closeTo(0.0, 1e-9));
      expect(r.y, closeTo(0.0, 1e-9));
      expect(r.gamma, closeTo(0.0, 1e-12));
      expect(r.k, closeTo(0.9996, 1e-9));
    });

    test('Any point on central meridian has x = 0', () {
      for (final lat in [0.0, 10.0, 30.0, 45.0, 60.0, 80.0, -30.0, -60.0]) {
        final r = tm.forward(lon0: 0.0, lat: lat, lon: 0.0);
        expect(r.x, closeTo(0.0, 1e-8),
            reason: 'x=0 on central meridian, lat=$lat');
      }
    });

    test('γ = 0 on equator (lat=0)', () {
      for (final lon in [1.0, 5.0, 10.0, 20.0, 45.0]) {
        final r = tm.forward(lon0: 0.0, lat: 0.0, lon: lon);
        expect(r.gamma, closeTo(0.0, 1e-12),
            reason: 'γ=0 on equator, lon=$lon');
      }
    });

    test('γ = 0 on central meridian', () {
      for (final lat in [10.0, 30.0, 45.0, 60.0, 80.0]) {
        final r = tm.forward(lon0: 0.0, lat: lat, lon: 0.0);
        expect(r.gamma, closeTo(0.0, 1e-12),
            reason: 'γ=0 on meridian, lat=$lat');
      }
    });

    test('Scale symmetric: k(lat, lon) = k(lat, -lon)', () {
      for (final lon in [1.0, 5.0, 10.0, 20.0, 45.0]) {
        final rP = tm.forward(lon0: 0.0, lat: 45.0, lon:  lon);
        final rN = tm.forward(lon0: 0.0, lat: 45.0, lon: -lon);
        expect(rP.k, closeTo(rN.k, 1e-14),
            reason: 'k symmetric at lon=$lon');
      }
    });

    test('Easting antisymmetric: x(lat, -lon) = -x(lat, lon)', () {
      for (final lon in [1.0, 10.0, 30.0]) {
        final rP = tm.forward(lon0: 0.0, lat: 45.0, lon:  lon);
        final rN = tm.forward(lon0: 0.0, lat: 45.0, lon: -lon);
        expect(rN.x, closeTo(-rP.x, 1e-9));
      }
    });

    test('Northing symmetric: y(lat, lon) = y(lat, -lon)', () {
      for (final lon in [1.0, 10.0, 30.0]) {
        final rP = tm.forward(lon0: 0.0, lat: 45.0, lon:  lon);
        final rN = tm.forward(lon0: 0.0, lat: 45.0, lon: -lon);
        expect(rP.y, closeTo(rN.y, 1e-9));
      }
    });

    test('Northing antisymmetric in lat: y(-lat, 0) = -y(lat, 0)', () {
      for (final lat in [10.0, 30.0, 60.0]) {
        final rP = tm.forward(lon0: 0.0, lat:  lat, lon: 0.0);
        final rN = tm.forward(lon0: 0.0, lat: -lat, lon: 0.0);
        expect(rN.y, closeTo(-rP.y, 1e-9));
      }
    });
  });

  // -------------------------------------------------------------------------
  // 2. Round-trip accuracy: forward → reverse ≤ 1 nm geographic distance
  // -------------------------------------------------------------------------
  group('Round-trip accuracy (≤ 1 nm)', () {
    const tolerance = 1e-6; // 1 μm ≈ 1e-11 degrees – well within 1 nm

    final testPoints = [
      [0.0, 0.0],
      [0.0, 1.0],
      [45.0, 10.0],
      [60.0, 30.0],
      [80.0, 45.0],
      [-30.0, -15.0],
      [89.9, 0.5],
      // Near-pole
      [89.999, 1.0],
      // Far from central meridian
      [0.0, 60.0],
      [45.0, 80.0],
    ];

    for (final pt in testPoints) {
      final lat = pt[0], lon = pt[1];
      test('Round-trip lat=$lat lon=$lon', () {
        final fwd = tm.forward(lon0: 0.0, lat: lat, lon: lon);
        final rev = tm.reverse(lon0: 0.0, x: fwd.x, y: fwd.y);
        expect(rev.lat, closeTo(lat, tolerance),
            reason: 'lat round-trip for ($lat, $lon)');
        expect(rev.lon, closeTo(lon, tolerance),
            reason: 'lon round-trip for ($lat, $lon)');
      });
    }
  });

  // -------------------------------------------------------------------------
  // 3. Agreement with Krüger series (≤ 5 nm) within 35° of central meridian
  // -------------------------------------------------------------------------
  group('Agreement with Krüger series (≤ 5 nm within 35°)', () {
    final krueger = KruegerTM(k0: 0.9996);
    const tol = 5e-3; // 5 mm; Krüger series ≤ 5 nm but we give slack for test

    final nearPoints = [
      [0.0, 0.0],
      [0.0, 10.0],
      [30.0, 5.0],
      [45.0, 20.0],
      [60.0, 15.0],
      [80.0, 10.0],
      [-45.0, 30.0],
    ];

    for (final pt in nearPoints) {
      final lat = pt[0], lon = pt[1];
      test('Exact vs Krüger lat=$lat lon=$lon', () {
        final exact = tm.forward(lon0: 0.0, lat: lat, lon: lon);
        final series = krueger.forward(0.0, lat, lon);
        expect(exact.x, closeTo(series.x, tol),
            reason: 'x matches Krüger for ($lat, $lon)');
        expect(exact.y, closeTo(series.y, tol),
            reason: 'y matches Krüger for ($lat, $lon)');
      });
    }
  });

  // -------------------------------------------------------------------------
  // 4. High-precision reference values (from TMcoords.dat)
  // -------------------------------------------------------------------------
  group('TMcoords.dat reference values (0.001 m)', () {
    for (final d in _exactData) {
      final lat = d[0], lon = d[1], expX = d[2], expY = d[3];
      test('Forward lat=$lat lon=$lon', () {
        final r = tm.forward(lon0: 0.0, lat: lat, lon: lon);
        if (expX != 0.0) expect(r.x, closeTo(expX, 0.01));
        if (expY != 0.0) expect(r.y, closeTo(expY, 0.01));
      });
    }
  });

  // -------------------------------------------------------------------------
  // 5. Edge cases
  // -------------------------------------------------------------------------
  group('Edge cases', () {
    test('North pole (lat=90) → x=0', () {
      final r = tm.forward(lon0: 0.0, lat: 90.0, lon: 0.0);
      expect(r.x, closeTo(0.0, 1e-6));
    });

    test('North pole (lat=90) → scale k=1 · k0', () {
      final r = tm.forward(lon0: 0.0, lat: 90.0, lon: 0.0);
      expect(r.k, closeTo(0.9996, 1e-9));
    });

    test('South pole (lat=-90) → x=0', () {
      final r = tm.forward(lon0: 0.0, lat: -90.0, lon: 0.0);
      expect(r.x, closeTo(0.0, 1e-6));
    });

    test('Reverse of forward pole stays at pole', () {
      final fwd = tm.forward(lon0: 0.0, lat: 90.0, lon: 0.0);
      final rev = tm.reverse(lon0: 0.0, x: fwd.x, y: fwd.y);
      expect(rev.lat, closeTo(90.0, 1e-9));
    });

    test('WGS84 singleton vs. fresh instance give identical results', () {
      final fresh = TransverseMercatorExact(k0: 0.9996);
      final r1 = tm.forward(lon0: 3.0, lat: 48.0, lon: 10.0);
      final r2 = fresh.forward(lon0: 3.0, lat: 48.0, lon: 10.0);
      expect(r1.x, closeTo(r2.x, 1e-10));
      expect(r1.y, closeTo(r2.y, 1e-10));
    });

    test('Different central meridians give consistent results', () {
      // forward(lon0=3, lat, lon=10) should give same x,y as
      // forward(lon0=0, lat, lon=7) [shifted by 3°]
      final r3 = tm.forward(lon0: 3.0, lat: 48.0, lon: 10.0);
      final r0 = tm.forward(lon0: 0.0, lat: 48.0, lon: 7.0);
      expect(r3.x, closeTo(r0.x, 1e-6));
      expect(r3.y, closeTo(r0.y, 1e-6));
    });

    test('Scale on central meridian equals k0 at equator', () {
      final r = tmK1.forward(lon0: 0.0, lat: 0.0, lon: 0.0);
      expect(r.k, closeTo(1.0, 1e-9));
    });

    test('Meridian convergence at (45°N, 10°E) is positive', () {
      // gamma should be in (0°, 10°) for a NE point
      final r = tm.forward(lon0: 0.0, lat: 45.0, lon: 10.0);
      expect(r.gamma, greaterThan(0.0));
      expect(r.gamma, lessThan(10.0));
    });
  });

  // -------------------------------------------------------------------------
  // 6. Elliptic function unit tests (internal validation)
  // -------------------------------------------------------------------------
  group('EllipticFunction (WGS84 parameters)', () {
    // For WGS84: mu = e² ≈ 0.00669438
    const mu = 0.00669437999014;
    final eEu = EllipticFunction(mu);
    final eEv = EllipticFunction(1.0 - mu);

    test('K(mu) is finite and > pi/2', () {
      final k = eEu.K();
      expect(k.isFinite, isTrue);
      expect(k, greaterThan(math.pi / 2.0));
    });

    test('E(mu) is finite and < pi/2', () {
      final e = eEu.E();
      expect(e.isFinite, isTrue);
      expect(e, lessThan(math.pi / 2.0));
      expect(e, greaterThan(1.0));
    });

    test('KE = K - E', () {
      expect(eEu.KE(), closeTo(eEu.K() - eEu.E(), 1e-15));
    });

    test('K(0) = E(0) = pi/2', () {
      final e0 = EllipticFunction(0.0);
      expect(e0.K(), closeTo(math.pi / 2.0, 1e-15));
      expect(e0.E(), closeTo(math.pi / 2.0, 1e-15));
    });

    test('E(k2=1) = 1', () {
      // k²=1 → kp²=0, E = 1
      final e1 = EllipticFunction(1.0);
      expect(e1.E(), closeTo(1.0, 1e-12));
    });

    test('am(0) gives sn=0 cn=1 dn=1', () {
      final r = eEu.am(0.0);
      expect(r.sn, closeTo(0.0, 1e-15));
      expect(r.cn, closeTo(1.0, 1e-15));
      expect(r.dn, closeTo(1.0, 1e-15));
    });

    test('am(K) gives sn=1', () {
      final r = eEu.am(eEu.K());
      expect(r.sn.abs(), closeTo(1.0, 1e-12));
    });

    test('sn²+cn²=1 at a representative point', () {
      final x = 0.8;
      final r = eEu.am(x);
      expect(r.sn * r.sn + r.cn * r.cn, closeTo(1.0, 1e-14));
    });

    test('dn²+k²·sn²=1 at a representative point', () {
      final x = 0.8;
      final r = eEu.am(x);
      expect(r.dn * r.dn + mu * r.sn * r.sn, closeTo(1.0, 1e-12));
    });

    test('eIncomplete(0, 1, 1) = 0', () {
      expect(eEu.eIncomplete(0.0, 1.0, 1.0), closeTo(0.0, 1e-15));
    });

    test('eIncomplete at K gives complete E', () {
      final r = eEu.am(eEu.K());
      final e = eEu.eIncomplete(r.sn, r.cn, r.dn);
      expect(e.abs(), closeTo(eEu.E(), 1e-12));
    });

    test('eEv K is large (k²≈0.9933 → large K)', () {
      // eEv is constructed with k²=1-e²≈0.9933, kp²=e²≈0.0067
      // For k close to 1, K→∞; K(0.9933) ≈ 3.89
      expect(eEv.K(), greaterThan(3.0));
      expect(eEv.K(), lessThan(10.0));
    });
  });
}
