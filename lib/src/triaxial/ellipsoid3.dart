// ellipsoid3.dart
// Dart port of GeographicLib::Triaxial::Ellipsoid3
// from Ellipsoid3.{hpp,cpp} / Ellipsoid3.js
//
// ==========================================================================
// Original C++/JS: Copyright (c) Charles Karney (2024-2025)
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
import '../geo_math.dart';
import 'angle.dart';

// ignore_for_file: non_constant_identifier_names

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

double _sq(double x) => x * x;

/// Normalize a 3-vector in place. Returns the length.
double _normvec(List<double> R) {
  final h = math.sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
  R[0] /= h;
  R[1] /= h;
  R[2] /= h;
  return h;
}

// ---------------------------------------------------------------------------
// Ellipsoid3 class
// ---------------------------------------------------------------------------

/// A triaxial ellipsoid with semiaxes a ≥ b ≥ c > 0.
class Ellipsoid3 {
  late final double _a;
  late final double _b;
  late final double _c;
  late final double _e2;
  late final double _k2;
  late final double _kp2;
  late final double _k;
  late final double _kp;
  late final bool _oblate;
  late final bool _prolate;
  late final bool _biaxial;

  // -------------------------------------------------------------------------
  // Constructors
  // -------------------------------------------------------------------------

  /// Construct from semiaxes [a] ≥ [b] ≥ [c] > 0.
  Ellipsoid3(double a, double b, double c) {
    _a = a;
    _b = b;
    _c = c;
    final s = (a - c) * (a + c);
    _e2 = s / _sq(b);
    if (s == 0.0) {
      _kp2 = 0.0;
      _k2 = 1.0;
    } else {
      _kp2 = (a - b) * (a + b) / s;
      _k2 = (b - c) * (b + c) / s;
    }
    _k = math.sqrt(_k2);
    _kp = math.sqrt(_kp2);
    _oblate = _kp2 == 0.0;
    _prolate = _k2 == 0.0;
    _biaxial = _oblate || _prolate;
  }

  /// Construct from shape parameters (b, e2, k2, kp2).
  ///
  /// [b]   – intermediate semiaxis
  /// [e2]  – eccentricity squared
  /// [k2]  – shape parameter k²
  /// [kp2] – complementary shape parameter k'²
  Ellipsoid3.fromParams(double b, double e2, double k2, double kp2) {
    _b = b;
    _e2 = e2;
    final ksum = k2 + kp2;
    _k2 = k2 / ksum;
    _kp2 = kp2 / ksum;
    _k = math.sqrt(_k2);
    _kp = math.sqrt(_kp2);
    _a = b * math.sqrt(1.0 + e2 * _kp2);
    _c = b * math.sqrt(1.0 - e2 * _k2);
    _oblate = _kp2 == 0.0;
    _prolate = _k2 == 0.0;
    _biaxial = _oblate || _prolate;
  }

  // -------------------------------------------------------------------------
  // Inspectors
  // -------------------------------------------------------------------------

  double get a => _a;
  double get b => _b;
  double get c => _c;
  double get e2 => _e2;
  double get k2 => _k2;
  double get kp2 => _kp2;
  double get k => _k;
  double get kp => _kp;
  bool get oblate => _oblate;
  bool get prolate => _prolate;
  bool get biaxial => _biaxial;

  // -------------------------------------------------------------------------
  // Static: Flip(bet, omg, alp)
  // bet.reflect(flipc)  → negate c only
  // omg.reflect(flips)  → negate s only
  // alp.reflect(flips, flipc) → negate both s and c
  // -------------------------------------------------------------------------

  /// Apply a coordinate flip to [bet], [omg], and optionally [alp].
  static void flip(Angle bet, Angle omg, [Angle? alp]) {
    bet.reflect(false, true);  // negate c of bet
    omg.reflect(true, false);  // negate s of omg
    alp?.reflect(true, true);  // negate both s and c of alp
  }

  // -------------------------------------------------------------------------
  // Static: AngNorm
  // Normalize so that bet in [-π/2, π/2] (alt=false)
  //                 or omg.s in [0, π] (alt=true).
  // Returns whether a flip was applied.
  // -------------------------------------------------------------------------

  /// Normalize coordinates.
  ///
  /// Two overloaded forms:
  ///   angNorm(bet, omg, alp, alt)  — with explicit alp Angle
  ///   angNorm(bet, omg, alt)       — without alp (pass alp=null)
  static bool angNorm(Angle bet, Angle omg, [Angle? alp, bool alt = false]) {
    bool flipF;
    if (alt) {
      flipF = bet.cx < 0; // 1/omg.s < 0 equivalent: sign of omg.sx negative
      // Actually: alt=true means normalize so omg.s >= 0
      // JS: flip = alt ? (1/omg.s() < 0) : (1/bet.c() < 0)
      // We need the signbit of omg.sx:
      flipF = omg.sx.isNegative || (omg.sx == 0.0 && (1.0 / omg.sx).isNegative);
    } else {
      flipF = bet.cx.isNegative || (bet.cx == 0.0 && (1.0 / bet.cx).isNegative);
    }
    if (flipF) {
      flip(bet, omg, alp);
    }
    return flipF;
  }

  // -------------------------------------------------------------------------
  // Norm(R, [V])
  // Scale R to surface of ellipsoid, optionally project V onto tangent plane.
  // -------------------------------------------------------------------------

  /// Scale Cartesian vector [R] to the ellipsoid surface.
  /// Optionally projects velocity [V] onto the tangent plane.
  void norm(List<double> R, [List<double>? V]) {
    final a = _a, b = _b, c = _c;
    final ra = math.sqrt(
        _sq(R[0] / a) + _sq(R[1] / b) + _sq(R[2] / c));
    R[0] /= ra; R[1] /= ra; R[2] /= ra;
    if (V != null) {
      final up = [R[0] / _sq(a), R[1] / _sq(b), R[2] / _sq(c)];
      final u2 = _sq(up[0]) + _sq(up[1]) + _sq(up[2]);
      final uv = up[0] * V[0] + up[1] * V[1] + up[2] * V[2];
      final f = uv / u2;
      V[0] -= f * up[0]; V[1] -= f * up[1]; V[2] -= f * up[2];
      _normvec(V);
    }
  }

  // -------------------------------------------------------------------------
  // cart2toellipint — internal helper
  // -------------------------------------------------------------------------

  ({Angle bet, Angle omg}) _cart2ToEllipInt(
      List<double> R, List<double> axes) {
    final a = axes[0], bAx = axes[1], c = axes[2];
    final k2 = _k2, kp2 = _kp2, k = _k, kp = _kp;
    final xi = R[0] / a, eta = R[1] / bAx, zeta = R[2] / c;
    var g = k2 * _sq(xi) + (k2 - kp2) * _sq(eta) - kp2 * _sq(zeta);
    if (R[0].abs() == a * kp2 && R[1] == 0.0 && R[2].abs() == c * k2) {
      g = 0.0;
    }
    final h = math.sqrt(g * g + _sq(2.0 * k * kp * eta));
    double so, co, sb, cb;
    if (h == 0.0) {
      so = 0.0; cb = 0.0;
    } else if (g < 0) {
      so = GeoMath.copysign(math.sqrt((h - g) / 2.0) / kp, eta);
      cb = (eta / so).abs();
    } else {
      cb = math.sqrt((h + g) / 2.0) / k;
      so = eta / cb;
    }
    final tz = math.sqrt(k * k + _sq(kp * so));
    final tx = math.sqrt(_sq(k * cb) + kp * kp);
    sb = tz == 0.0 ? -1.0 : zeta / tz;
    co = tx == 0.0 ? 1.0 : xi / tx;
    final bet = Angle(sb, cb, 0, true);
    final omg = Angle(so, co, 0, true);
    return (bet: bet, omg: omg);
  }

  // -------------------------------------------------------------------------
  // Public: cart2toellip(R) → {bet, omg}
  // -------------------------------------------------------------------------

  /// Convert Cartesian point [R] to elliptic coordinates (bet, omg).
  ({Angle bet, Angle omg}) cart2toellip(List<double> R) {
    return _cart2ToEllipInt(R, [_a, _b, _c]);
  }

  // -------------------------------------------------------------------------
  // cart2toellipDir(bet, omg, V) → alp
  // Convert cartesian direction V to heading alp.
  // -------------------------------------------------------------------------

  /// Convert Cartesian direction [V] to heading angle [alp].
  Angle cart2toellipDir(Angle bet, Angle omg, List<double> V) {
    final a = _a, b = _b, c = _c;
    final k2 = _k2, kp2 = _kp2;
    final tz = math.sqrt(_sq(_k) + _sq(_kp * omg.s()));
    final tx = math.sqrt(_sq(_k * bet.c()) + _sq(_kp));
    late Angle alp;
    if (tx == 0.0 || tz == 0.0 || !(bet.c() == 0.0 && omg.s() == 0.0)) {
      // Not a triaxial umbilical point
      List<double> N, E;
      if (tx == 0.0) {
        final scb = bet.cx < 0 ? -1.0 : 1.0;
        N = [-omg.c() * bet.s() * scb, -omg.s() * bet.s(), tx * bet.s()];
        E = [-omg.s(), omg.c() * scb, tx];
      } else if (tz == 0.0) {
        final sso = omg.sx.isNegative ? -1.0 : 1.0;
        N = [tz, -bet.s() * sso, bet.c()];
        E = [tz * omg.c(), bet.c() * omg.c(), bet.s() * omg.c() * sso];
      } else {
        N = [
          -a * k2 * bet.c() * bet.s() * omg.c() / tx,
          -b * bet.s() * omg.s(),
          c * bet.c() * tz
        ];
        E = [
          -a * tx * omg.s(),
          b * bet.c() * omg.c(),
          c * kp2 * bet.s() * omg.c() * omg.s() / tz
        ];
      }
      _normvec(N); _normvec(E);
      alp = Angle(
          V[0] * E[0] + V[1] * E[1] + V[2] * E[2],
          V[0] * N[0] + V[1] * N[1] + V[2] * N[2],
          0, true);
    } else {
      // Special treatment at umbilical points
      final w = bet.s() * omg.c();
      var upx = omg.c() * tx / a;
      var upz = bet.s() * tz / c;
      final h2 = math.sqrt(upx * upx + upz * upz);
      if (h2 != 0.0) { upx /= h2; upz /= h2; }
      final s2a = -V[1] * w;
      final c2a = (upz * V[0] - upx * V[2]) * w;
      final flipV = -bet.s();
      if (c2a >= 0) {
        alp = Angle(flipV * s2a, flipV * (1.0 + c2a), 0, false);
      } else {
        alp = Angle(
            flipV * GeoMath.copysign(1.0 - c2a, s2a),
            flipV * s2a.abs(),
            0, false);
      }
    }
    return alp;
  }

  // -------------------------------------------------------------------------
  // elliptocart2(bet, omg) → R
  // -------------------------------------------------------------------------

  /// Convert elliptic coordinates (bet, omg) to Cartesian [R].
  List<double> elliptocart2(Angle bet, Angle omg) {
    final a = _a, b = _b, c = _c;
    final k = _k, kp = _kp;
    final tx = math.sqrt(_sq(k * bet.c()) + _sq(kp));
    final tz = math.sqrt(_sq(k) + _sq(kp * omg.s()));
    return [a * omg.c() * tx, b * bet.c() * omg.s(), c * bet.s() * tz];
  }

  // -------------------------------------------------------------------------
  // elliptocart2dir(bet, omg, alp) → {R, V}
  // -------------------------------------------------------------------------

  /// Convert elliptic coordinates + heading to Cartesian position [R] and
  /// direction [V].
  ({List<double> R, List<double> V}) elliptocart2dir(
      Angle bet, Angle omg, Angle alp) {
    final a = _a, b = _b, c = _c;
    final k = _k, kp = _kp, k2 = _k2, kp2 = _kp2;
    final tx = math.sqrt(_sq(k * bet.c()) + _sq(kp));
    final tz = math.sqrt(_sq(k) + _sq(kp * omg.s()));
    final R = [a * omg.c() * tx, b * bet.c() * omg.s(), c * bet.s() * tz];
    List<double> V;
    if (bet.c() == 0.0 && omg.s() == 0.0 && !(k == 0.0 || kp == 0.0)) {
      // umbilical point
      final sa2 = 2.0 * alp.s() * alp.c();
      final ca2 = (alp.c() - alp.s()) * (alp.c() + alp.s());
      V = [
        a * k / b * omg.c() * ca2,
        -omg.c() * bet.s() * sa2,
        -c * kp / b * bet.s() * ca2
      ];
    } else {
      List<double> N, E;
      if (tx == 0.0) {
        final scb = bet.cx < 0 ? -1.0 : 1.0;
        N = [-omg.c() * bet.s() * scb, -omg.s() * bet.s(), 0.0];
        E = [-omg.s(), omg.c() * scb, 0.0];
      } else if (tz == 0.0) {
        final sso = omg.sx.isNegative ? -1.0 : 1.0;
        N = [0.0, -bet.s() * sso, bet.c()];
        E = [0.0, bet.c() * omg.c(), bet.s() * omg.c() * sso];
      } else {
        N = [
          -a * k2 * bet.c() * bet.s() * omg.c() / tx,
          -b * bet.s() * omg.s(),
          c * bet.c() * tz
        ];
        E = [
          -a * tx * omg.s(),
          b * bet.c() * omg.c(),
          c * kp2 * bet.s() * omg.c() * omg.s() / tz
        ];
      }
      _normvec(N); _normvec(E);
      V = [
        alp.c() * N[0] + alp.s() * E[0],
        alp.c() * N[1] + alp.s() * E[1],
        alp.c() * N[2] + alp.s() * E[2]
      ];
    }
    return (R: R, V: V);
  }
}

