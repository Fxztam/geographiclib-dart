// geo_math.dart
// Dart port of Math.js from geographiclib-geodesic v2.2.0
//
// ==========================================================================
// Original C++/JS: Copyright (c) Charles Karney (2011-2021)
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

/// Utility math functions for geodesic computations.
class GeoMath {
  /// Number of mantissa bits in a double.
  static const int digits = 53;

  /// Machine epsilon (2^-52).
  static const double epsilon = 2.220446049250313e-16; // pow(0.5, 52)

  /// Factor to convert degrees to radians.
  static const double degree = math.pi / 180.0;

  /// Square.
  static double sq(double x) => x * x;

  /// Hypotenuse (avoids Math.hypot inaccuracy on some platforms).
  static double hypot(double x, double y) => math.sqrt(x * x + y * y);

  /// Real cube root.
  static double cbrt(double x) {
    final y = math.pow(x.abs(), 1.0 / 3.0) as double;
    return x > 0 ? y : (x < 0 ? -y : x);
  }

  /// log(1 + x) accurately.
  static double log1p(double x) {
    final y = 1.0 + x;
    final z = y - 1.0;
    return z == 0.0 ? x : x * math.log(y) / z;
  }

  /// Inverse hyperbolic tangent.
  static double atanh(double x) {
    final y = x.abs();
    final r = log1p(2.0 * y / (1.0 - y)) / 2.0;
    return x > 0 ? r : (x < 0 ? -r : x);
  }

  /// Copy the sign of [y] onto the magnitude of [x].
  static double copysign(double x, double y) {
    return x.abs() * ((y < 0 || (y == 0.0 && (1.0 / y).isNegative)) ? -1.0 : 1.0);
  }

  /// Error-free sum. Returns [s, t] where s = round(u+v) and t = u+v - s.
  static (double s, double t) sum(double u, double v) {
    final s = u + v;
    final up = s - v;
    final vpp = s - up;
    final t = s != 0.0 ? 0.0 - ((up - u) + (vpp - v)) : s;
    return (s, t);
  }

  /// Evaluate polynomial using Horner's method.
  /// [N] is degree, [p] is coefficient array (leading first), [s] is start index, [x] is variable.
  static double polyval(int N, List<double> p, int s, double x) {
    double y = N < 0 ? 0.0 : p[s++];
    while (--N >= 0) {
      y = y * x + p[s++];
    }
    return y;
  }

  /// Coarsen a value close to zero to avoid near-singular cases.
  static double angRound(double x) {
    const z = 1.0 / 16.0;
    if (x.isNaN) return x;
    final y = x.abs();
    // The compiler mustn't simplify z - (z - y) to y
    final yr = y < z ? z - (z - y) : y;
    return copysign(yr, x);
  }

  /// Remainder in range [-y/2, y/2] — IEEE 754 semantics.
  /// Rounds to nearest even on ties; result is 0.0 with sign of x.
  static double remainder(double x, double y) {
    if (x.isNaN || x.isInfinite || y == 0.0 || y.isNaN) return double.nan;
    // Compute n = round(x/y) using round-half-to-even (banker's rounding).
    final q = x / y;
    final t = q.truncateToDouble();
    final frac = (q - t).abs();
    final double n;
    if (frac != 0.5) {
      n = q.roundToDouble();
    } else {
      // Tie: choose nearest even integer.
      n = t.remainder(2.0) == 0.0 ? t : t + (q >= 0 ? 1.0 : -1.0);
    }
    final r = x - n * y;
    // IEEE 754: if result is zero, its sign equals sign of x.
    return r == 0.0 ? copysign(0.0, x) : r;
  }

  /// Normalize angle to [-180, 180], preserving sign of ±0 and ±180.
  static double angNormalize(double x) {
    final y = remainder(x, 360.0);
    return y.abs() == 180.0 ? copysign(180.0, x) : y;
  }

  /// Return NaN if |x| > 90, else return x.
  static double latFix(double x) => x.abs() > 90.0 ? double.nan : x;

  /// Exact difference of two angles: returns (d, e) where d+e = y-x mod 360°,
  /// d is the nearest representable float, and e is the rounding error.
  static (double d, double e) angDiff(double x, double y) {
    final r1 = sum(remainder(-x, 360.0), remainder(y, 360.0));
    final r2 = sum(remainder(r1.$1, 360.0), r1.$2);
    double d = r2.$1;
    final double e = r2.$2;
    if (d == 0.0 || d.abs() == 180.0) {
      d = copysign(d, e == 0.0 ? y - x : -e);
    }
    return (d, e);
  }

  /// sin and cos of angle in degrees with exact results at multiples of 30° and 45°.
  /// Signed-zero: sincosd(-0.0) has sin = -0.0.
  static (double s, double c) sincosd(double x) {
    if (x.isNaN || x.isInfinite) return (double.nan, double.nan);
    double d = x.remainder(360.0); // truncate-toward-zero (matches JS % semantics)
    final q = (d / 90.0).roundToDouble().toInt();
    d -= 90.0 * q;
    // now abs(d) <= 45
    final r = d * degree;
    double s = math.sin(r);
    double c = math.cos(r);
    if (d.abs() == 45.0) {
      c = math.sqrt(0.5);
      s = copysign(c, r);
    } else if (d.abs() == 30.0) {
      c = math.sqrt(0.75);
      s = copysign(0.5, r);
    }
    double sinx, cosx;
    switch (q & 3) {
      case 0:  sinx =  s; cosx =  c; break;
      case 1:  sinx =  c; cosx = -s; break;
      case 2:  sinx = -s; cosx = -c; break;
      default: sinx = -c; cosx =  s; break;
    }
    cosx = cosx + 0.0; // ensure +0.0
    if (sinx == 0.0) sinx = copysign(sinx, x);
    return (sinx, cosx);
  }

  /// sin and cos of (x + t) degrees where x is already reduced and t is a small correction.
  static (double s, double c) sincosde(double x, double t) {
    if (x.isNaN || x.isInfinite) return (double.nan, double.nan);
    double d = x.remainder(360.0); // truncate-toward-zero (matches JS % semantics)
    final q = (d / 90.0).roundToDouble().toInt();
    d = angRound((d - 90.0 * q) + t);
    final r = d * degree;
    double s = math.sin(r);
    double c = math.cos(r);
    if (d.abs() == 45.0) {
      c = math.sqrt(0.5);
      s = copysign(c, r);
    } else if (d.abs() == 30.0) {
      c = math.sqrt(0.75);
      s = copysign(0.5, r);
    }
    double sinx, cosx;
    switch (q & 3) {
      case 0:  sinx =  s; cosx =  c; break;
      case 1:  sinx =  c; cosx = -s; break;
      case 2:  sinx = -s; cosx = -c; break;
      default: sinx = -c; cosx =  s; break;
    }
    cosx = cosx + 0.0;
    if (sinx == 0.0) sinx = copysign(sinx, x + t);
    return (sinx, cosx);
  }

  /// atan2 in degrees with exact results for ±∞ and signed-zero handling.
  static double atan2d(double y, double x) {
    int q = 0;
    double yt = y, xt = x;
    if (yt.abs() > xt.abs()) {
      final tmp = xt; xt = yt; yt = tmp;
      q = 2;
    }
    if (copysign(1.0, xt) < 0) { xt = -xt; ++q; }
    double ang = math.atan2(yt, xt) / degree;
    switch (q) {
      case 1: ang = copysign(180.0, y) - ang; break;
      case 2: ang = 90.0 - ang; break;
      case 3: ang = -90.0 + ang; break;
      default: break;
    }
    return ang;
  }
}

/// Compensated accumulator for summing many numbers with double precision.
class Accumulator {
  double _s;
  double _t;

  Accumulator([dynamic y]) : _s = 0.0, _t = 0.0 {
    if (y == null) {
      _s = 0.0; _t = 0.0;
    } else if (y is Accumulator) {
      _s = y._s; _t = y._t;
    } else {
      _s = (y as double) + 0.0; _t = 0.0;
    }
  }

  void set(dynamic y) {
    if (y is Accumulator) {
      _s = y._s; _t = y._t;
    } else {
      _s = (y as double) + 0.0; _t = 0.0;
    }
  }

  void add(double y) {
    final u = GeoMath.sum(y, _t);
    final v = GeoMath.sum(u.$1, _s);
    final ut = u.$2;
    _s = v.$1;
    _t = v.$2;
    if (_s == 0.0) {
      _s = ut;
    } else {
      _t += ut;
    }
  }

  /// Return accumulated sum optionally plus [y].
  double sum([double? y]) {
    if (y == null) return _s;
    final b = Accumulator(this);
    b.add(y);
    return b._s;
  }

  void negate() {
    _s *= -1.0;
    _t *= -1.0;
  }

  void remainder(double y) {
    _s = GeoMath.remainder(_s, y);
    add(0.0);
  }
}
