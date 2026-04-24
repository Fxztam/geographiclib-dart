// angle.dart
// Dart port of GeographicLib::AngleT<double> from Angle.hpp / Angle.js
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

// ignore_for_file: non_constant_identifier_names

// ---------------------------------------------------------------------------
// Module-level constants
// ---------------------------------------------------------------------------

const double _pi = math.pi;
const double _td = 360.0; // full period in degrees
const double _eps = GeoMath.epsilon; // machine epsilon ≈ 2^-52
const double _z = 1.0 / 1024.0; // for rnd()
const double _epsAngle = _eps / (1 << 20); // tiny angle

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Round to nearest integer (same as C++ std::rint).
/// Returns 0 for non-finite inputs (guards against NaN/Infinity turn counts).
int _rint(double x) {
  if (!x.isFinite) return 0;
  return x.roundToDouble().toInt();
}

/// Test whether x is negative-zero or negative finite.
/// Equivalent to JS `1/x < 0` (signbit).
bool _signbit(double x) =>
    x.isNegative || (x == 0.0 && (1.0 / x).isNegative);

/// Flush tiny values to ±0 (matches AngleT::rnd).
/// z = 1/1024 so the smallest gap is eps/2048.
double _rnd(double x) {
  final y = x.abs();
  final w = _z - y;
  final yy = w > 0 ? _z - w : y;
  return GeoMath.copysign(yy, x);
}

// ---------------------------------------------------------------------------
// Angle class
// ---------------------------------------------------------------------------

/// Represents an angle as normalised (sin, cos) pair plus integer turn count.
///
/// Stores the angle θ as:
///   sx = sin(θ mod 2π)
///   cx = cos(θ mod 2π)
///   nx = number of complete turns
///
/// The pair (sx, cx) is always on the unit circle (or NaN for NaN angles).
class Angle {
  double sx;
  double cx;
  int nx;

  // Primary constructor
  // -------------------------------------------------------------------------

  /// Construct from sin, cos, and optional turn count.
  ///
  /// [s] – sine component
  /// [c] – cosine component
  /// [n] – number of complete turns (rounded to nearest integer)
  /// [normp] – if true, (s,c) is already normalised to the unit circle
  Angle(double s, double c, [int n = 0, bool normp = false])
      : sx = s,
        cx = c,
        nx = n {
    if (!normp) {
      final h = math.sqrt(sx * sx + cx * cx);
      if (h == 0.0) {
        cx = GeoMath.copysign(1.0, cx);
      } else if (h.isFinite) {
        sx /= h;
        cx /= h;
      } else if (h.isNaN || (!sx.isFinite && !cx.isFinite)) {
        sx = double.nan;
        cx = double.nan;
      } else if (!sx.isFinite) {
        sx = GeoMath.copysign(1.0, sx);
        cx = GeoMath.copysign(0.0, cx);
      } else {
        sx = GeoMath.copysign(0.0, sx);
        cx = GeoMath.copysign(1.0, cx);
      }
    }
  }

  // -------------------------------------------------------------------------
  // Default (zero angle): s=0, c=1, n=0
  // -------------------------------------------------------------------------

  /// The zero angle (0°).
  Angle.zero()
      : sx = 0.0,
        cx = 1.0,
        nx = 0;

  // -------------------------------------------------------------------------
  // Named factory constructors
  // -------------------------------------------------------------------------

  /// Create an Angle from [deg] degrees.
  factory Angle.fromDegrees(double deg) {
    final sc = GeoMath.sincosd(deg);
    final a = Angle(sc.$1, sc.$2, 0, true);
    a.nx = _rint((deg - GeoMath.atan2d(sc.$1, sc.$2)) / _td);
    return a;
  }

  /// Create an Angle from [rad] radians.
  factory Angle.fromRadians(double rad) {
    final sn = math.sin(rad);
    final cs = math.cos(rad);
    return Angle(sn, cs, _rint((rad - math.atan2(sn, cs)) / (2.0 * _pi)), true);
  }

  /// Create an Angle from a lambertian q = asinh(tan(angle)).
  factory Angle.fromLambertian(double q) {
    return Angle((math.exp(q) - math.exp(-q)) / 2.0, 1.0, 0);
  }

  /// Not-a-number angle.
  factory Angle.nan() {
    return Angle(double.nan, double.nan, 0, true);
  }

  /// Cardinal direction: [q] quarter-turns (rounded to nearest integer).
  /// +/-0 are distinguished.
  factory Angle.cardinal(int q) {
    if (!q.isFinite) return Angle.nan();
    final qi = q; // already int
    var iq = ((qi % 4) + 4) % 4; // positive modulo
    if (iq > 2) iq -= 4; // map to [-2, 2]
    double z = 0.0, s, c;
    switch (iq) {
      case -2:
        s = -z;
        c = -1.0;
        break;
      case -1:
        s = -1.0;
        c = z;
        break;
      case 1:
        s = 1.0;
        c = z;
        break;
      case 2:
        s = z;
        c = -1.0;
        break;
      default:
        // iq == 0: distinguish +0 vs -0
        s = (q != 0) ? z : 0.0;
        c = 1.0;
    }
    return Angle(s, c, (qi - iq) ~/ 4, true);
  }

  /// A tiny angle (epsilon/2^20 radians).
  factory Angle.epsilon() {
    return Angle(_epsAngle, 1.0, 0, true);
  }

  // -------------------------------------------------------------------------
  // Copy
  // -------------------------------------------------------------------------

  /// Return a copy of this angle.
  Angle copy() => Angle(sx, cx, nx, true);

  // -------------------------------------------------------------------------
  // Accessors
  // -------------------------------------------------------------------------

  /// Sine component.
  double s() => sx;

  /// Cosine component.
  double c() => cx;

  /// Tangent = s/c.
  double t() => sx / cx;

  /// Number of complete turns (converts -0 to +0).
  int n() => nx + 0;

  /// Turns adjusted for the -180° boundary case.
  int n0() {
    return (nx - (sx == 0.0 && _signbit(sx) && cx < 0 ? 1 : 0)) + 0;
  }

  // -------------------------------------------------------------------------
  // Conversions
  // -------------------------------------------------------------------------

  /// Convert to degrees (including turns).
  double degrees() {
    final d = degrees0();
    return nx == 0 ? d : d + _td * nx;
  }

  /// Convert to degrees, ignoring turns (result in (-180, 180]).
  double degrees0() => GeoMath.atan2d(sx, cx);

  /// Convert to radians (including turns).
  double radians() {
    final r = radians0();
    return nx == 0 ? r : r + 2.0 * _pi * nx;
  }

  /// Convert to radians, ignoring turns (result in (-π, π]).
  double radians0() => math.atan2(sx, cx);

  /// Return lambertian = asinh(tan(angle)).
  double lam() => math.log(sx / cx + math.sqrt(1 + (sx / cx) * (sx / cx)));
  // (math.log(x + sqrt(x^2+1)) = asinh(x))

  // -------------------------------------------------------------------------
  // Quadrant / cardinal helpers
  // -------------------------------------------------------------------------

  /// Nearest cardinal direction index (as a real number).
  double ncardinal() {
    final s = sx, c = cx;
    int iq;
    if (c < 0) {
      iq = (c <= -s.abs()) ? 2 : 1;
    } else {
      iq = (c >= s.abs()) ? 0 : 1;
    }
    iq = (_signbit(s) ? -1 : 1) * iq;
    return 4.0 * nx + iq;
  }

  /// Quadrant index 0–3.
  int quadrant() {
    final bs = _signbit(sx) ? 1 : 0;
    final bc = _signbit(cx) ? 1 : 0;
    return 2 * bs + (bc ^ bs);
  }

  /// Nearest cardinal direction as an Angle.
  Angle nearest([int ind = 0]) {
    double s, c;
    if (ind == 0) {
      if (cx.abs() >= sx.abs()) {
        s = GeoMath.copysign(0.0, sx);
        c = GeoMath.copysign(1.0, cx);
      } else {
        s = GeoMath.copysign(1.0, sx);
        c = GeoMath.copysign(0.0, cx);
      }
    } else if ((ind & 1) == 0) {
      s = GeoMath.copysign(0.0, sx);
      c = GeoMath.copysign(1.0, cx);
    } else {
      s = GeoMath.copysign(1.0, sx);
      c = GeoMath.copysign(0.0, cx);
    }
    return Angle(s, c, nx, true);
  }

  // -------------------------------------------------------------------------
  // Arithmetic operators (returning new Angle)
  // -------------------------------------------------------------------------

  /// Unary minus.
  Angle neg() => Angle(-sx, cx, -nx, true);

  /// Return this + [p].
  Angle add(Angle p) {
    final t = copy();
    return t..addEq(p);
  }

  /// Return this - [p].
  Angle sub(Angle p) {
    final t = copy();
    return t..subEq(p);
  }

  // -------------------------------------------------------------------------
  // In-place arithmetic
  // -------------------------------------------------------------------------

  /// Add [p] to this in place, keeping track of turns.
  void addEq(Angle p) {
    final q = ncardinal() + p.ncardinal();
    final c = cx * p.cx - sx * p.sx;
    sx = sx * p.cx + cx * p.sx;
    cx = c;
    nx += p.nx;
    final q2 = ncardinal();
    nx += _rint((q - q2) / 4.0);
  }

  /// Subtract [p] from this in place.
  void subEq(Angle p) => addEq(p.neg());

  // -------------------------------------------------------------------------
  // Comparison
  // -------------------------------------------------------------------------

  /// Test whether this ≈ 0 (within [mult]*epsilon).
  bool zerop([double mult = 0.0]) {
    return nx == 0 && cx > 0 && sx.abs() <= mult * _eps;
  }

  /// Equality: this == [p].
  bool eq(Angle p) => sub(p).zerop();

  // -------------------------------------------------------------------------
  // Modifying operations (all return this for chaining)
  // -------------------------------------------------------------------------

  /// Flush tiny s/c values to ±0.
  Angle round() {
    sx = _rnd(sx);
    cx = _rnd(cx);
    return this;
  }

  /// Re-normalise (s,c) to the unit circle.
  Angle renormalize() {
    final h = math.sqrt(sx * sx + cx * cx);
    sx /= h;
    cx /= h;
    return this;
  }

  /// Set the turn count to [n] (default 0).
  Angle setn([int n = 0]) {
    nx = n;
    return this;
  }

  /// Set turns so that the angle is in (-180°, +180°].
  Angle setn0([int n = 0]) {
    final extra = (sx == 0.0 && _signbit(sx) && cx < 0) ? 1 : 0;
    nx = n + extra;
    return this;
  }

  /// Set signs of s and c according to quadrant [q] (low 2 bits).
  /// q=0: (+,+), q=1: (+,-), q=2: (-,-), q=3: (-,+)
  Angle setquadrant(int q) {
    sx = GeoMath.copysign(sx, (q & 2) != 0 ? -1.0 : 1.0);
    cx = GeoMath.copysign(cx, ((q >> 1) ^ q) & 1 != 0 ? -1.0 : 1.0);
    return this;
  }

  /// Reflect / swap components.
  ///
  /// [flips] – negate s
  /// [flipc] – negate c
  /// [swapp] – swap s and c
  Angle reflect([bool flips = false, bool flipc = false, bool swapp = false]) {
    if (flips) sx *= -1.0;
    if (flipc) cx *= -1.0;
    if (swapp) {
      final tmp = sx;
      sx = cx;
      cx = tmp;
    }
    return this;
  }

  // -------------------------------------------------------------------------
  // Operations returning new Angle
  // -------------------------------------------------------------------------

  /// Return angle with turn count set to 0.
  Angle base() => Angle(sx, cx, 0, true);

  /// Return angle rebased to be within ±180° of center [cAng].
  Angle rebase(Angle cAng) {
    final t = copy();
    return t..setn0(sub(cAng).base().add(cAng).n0());
  }

  /// Return this if signbit(mult) is false, else return -this.
  Angle flipsign(double mult) {
    return _signbit(mult) ? neg() : copy();
  }

  /// Return atan(m * tan(this)), keeping turns continuous.
  /// [mv] must be ≥ 0.
  Angle modang(double mv) {
    if (_signbit(mv)) return Angle.nan(); // m < 0 → NaN
    return Angle(
      sx * (mv > 1.0 ? 1.0 : mv),
      cx / (mv > 1.0 ? mv : 1.0),
      nx,
    );
  }

  // -------------------------------------------------------------------------
  // Debug
  // -------------------------------------------------------------------------

  @override
  String toString() =>
      'Angle(s=$sx, c=$cx, n=$nx) = ${degrees().toStringAsFixed(6)}°';
}





