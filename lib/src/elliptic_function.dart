// elliptic_function.dart
// Dart port of GeographicLib::EllipticFunction (C++ → Dart)
//
// Provides Carlson symmetric-form elliptic integrals RF, RD, RG, RC, RJ,
// Jacobi elliptic functions sncndn and amplitude am, and the standard
// complete/incomplete integrals K, E, D, Pi, KE.
//
// ==========================================================================
// Original C++: Copyright (c) Charles Karney (2008-2024)
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
// Only the subset of the class required by TransverseMercatorExact is
// implemented:
//   EllipticFunction(k2)  - constructor with alpha2 = 0
//   K(), E(), KE()        – complete integrals
//   am(x, sn, cn, dn)     – Jacobi amplitude + sn/cn/dn via AGM
//   E(sn, cn, dn)         – incomplete E integral

import 'dart:math' as math;

/// Elliptic integrals and Jacobi elliptic functions.
///
/// Construct with a parameter [k2] = k² ∈ (−∞, 1].  The complementary
/// parameter is kp2 = 1 − k2 ≥ 0.
///
/// Only the methods required by [TransverseMercatorExact] are exposed:
/// [K], [E], [KE] (complete integrals) and [am], [eIncomplete] (incomplete).
class EllipticFunction {
  final double _k2;
  final double _kp2;

  late final double _kKc; // K(k)
  late final double _eEc; // E(k)

  // Tolerances derived from machine epsilon (double precision)
  static final double _tolRF =
      math.pow(3 * 2.220446049250313e-16 * 0.01, 1 / 8.0) as double;
  static final double _tolRD =
      math.pow(0.2 * (2.220446049250313e-16 * 0.01), 1 / 8.0) as double;
  static final double _tolAM =
      math.pow(2.220446049250313e-16, 0.75) as double;

  static const int _numIter = 13; // max AGM iterations (quadratic convergence)

  /// Creates an [EllipticFunction] for modulus-squared [k2] (alpha2 = 0).
  ///
  /// [k2] must satisfy k2 ≤ 1.  kp2 = 1 − k2 ≥ 0.
  EllipticFunction(double k2)
      : _k2 = k2,
        _kp2 = 1.0 - k2 {
    if (_k2 != 0.0) {
      _kKc = _kp2 != 0.0 ? _RF(_kp2, 1.0) : double.infinity;
      _eEc = _kp2 != 0.0 ? 2.0 * _RG(_kp2, 1.0) : 1.0;
    } else {
      _kKc = math.pi / 2.0;
      _eEc = math.pi / 2.0;
    }
  }

  // ---------------------------------------------------------------------------
  // Complete integrals
  // ---------------------------------------------------------------------------

  /// Complete elliptic integral K(k).
  double K() => _kKc;

  /// Complete elliptic integral E(k).
  double E() => _eEc;

  /// K(k) − E(k).
  double KE() => _kKc - _eEc;

  // ---------------------------------------------------------------------------
  // Jacobi amplitude and sn/cn/dn
  // ---------------------------------------------------------------------------

  /// Computes the Jacobi amplitude φ = am(x) and simultaneously stores
  /// sn = sin(φ), cn = cos(φ), dn = √(1 − k²·sin²(φ)) into the returned
  /// record.
  ///
  /// Uses the AGM algorithm (Sala 1989 / DLMF §22.20(ii)).
  ({double phi, double sn, double cn, double dn}) am(double x) {
    final phi = _am(x);
    double sn, cn, dn;
    if (_kp2 == 0.0) {
      sn = _tanh(x);
      cn = dn = 1.0 / _cosh(x);
    } else {
      sn = math.sin(phi);
      cn = math.cos(phi);
      dn = _delta(sn, cn);
    }
    return (phi: phi, sn: sn, cn: cn, dn: dn);
  }

  /// Incomplete elliptic integral E(sn, cn, dn).
  ///
  /// Carlson eq. 4.6 / DLMF 19.25.E9–E11.
  double eIncomplete(double sn, double cn, double dn) {
    final cn2 = cn * cn;
    final dn2 = dn * dn;
    final sn2 = sn * sn;
    double ei;
    if (cn2 != 0.0) {
      final absSn = sn.abs();
      if (_k2 <= 0.0) {
        // Carlson 4.6 / DLMF 19.25.E9
        ei = absSn * (_RF(cn2, dn2, 1.0) - _k2 * sn2 * _RD(cn2, dn2, 1.0) / 3.0);
      } else if (_kp2 >= 0.0) {
        // DLMF 19.25.E10
        ei = absSn * (_kp2 * _RF(cn2, dn2, 1.0) +
            _k2 * _kp2 * sn2 * _RD(cn2, 1.0, dn2) / 3.0 +
            _k2 * cn.abs() / dn);
      } else {
        // DLMF 19.25.E11
        ei = absSn * (-_kp2 * sn2 * _RD(dn2, 1.0, cn2) / 3.0 + dn / cn.abs());
      }
    } else {
      ei = E();
    }
    if (cn.isNegative) ei = 2.0 * E() - ei;
    return _copySign(ei, sn);
  }

  // ---------------------------------------------------------------------------
  // Private: Jacobi amplitude via AGM
  // ---------------------------------------------------------------------------

  double _am(double x) {
    if (_k2 == 0.0) return x;
    if (_kp2 == 0.0) return math.atan(_sinh(x)); // gd(x)

    double k2 = _k2, kp2 = _kp2;
    double xScaled = x;
    if (_k2 < 0.0) {
      // Sala Eq. 5.8
      k2 = -_k2 / _kp2; // used by the sqrt below
      kp2 = 1.0 / _kp2;
      xScaled *= math.sqrt(_kp2);
    }
    // suppress unused-variable lint for k2 computed above
    // ignore: unused_local_variable

    final a = List<double>.filled(_numIter, 0.0);
    final c = List<double>.filled(_numIter, 0.0);
    a[0] = 1.0;
    double b = math.sqrt(kp2);
    c[0] = math.sqrt(k2);
    int l = 1;
    for (; l < _numIter; ++l) {
      a[l] = (a[l - 1] + b) / 2.0;
      c[l] = (a[l - 1] - b) / 2.0;
      b = math.sqrt(a[l - 1] * b);
      if (!(c[l] > _tolAM * a[l])) break;
    }
    // a[l] ≈ π/(2K)
    double phi = a[l] * xScaled * (1 << l).toDouble();
    double phi1 = 0.0;
    for (; l > 0; --l) {
      phi1 = phi;
      phi = (phi + math.asin(c[l] * math.sin(phi) / a[l])) / 2.0;
    }
    return _k2 < 0.0 ? phi1 - phi : phi;
  }

  // Δ(sn, cn) = √(1 − k²·sn²)
  double _delta(double sn, double cn) =>
      math.sqrt(_kp2 + _k2 * cn * cn); // = sqrt(kp2 + k2*cn^2) = dn

  // ---------------------------------------------------------------------------
  // Carlson symmetric-form elliptic integrals (internal helpers)
  // ---------------------------------------------------------------------------

  /// RF(x, y, z) – Carlson eqs 2.2–2.7.
  static double _RF(double x, double y, [double? z]) {
    if (z == null) return _RF2(x, y); // degenerate overload
    double A0 = (x + y + z) / 3.0;
    double An = A0;
    double Q = math.max(math.max((A0 - x).abs(), (A0 - y).abs()), (A0 - z).abs()) / _tolRF;
    double x0 = x, y0 = y, z0 = z;
    double mul = 1.0;
    while (Q >= mul * An.abs()) {
      final lam = math.sqrt(x0) * math.sqrt(y0) +
          math.sqrt(y0) * math.sqrt(z0) +
          math.sqrt(z0) * math.sqrt(x0);
      An = (An + lam) / 4.0;
      x0 = (x0 + lam) / 4.0;
      y0 = (y0 + lam) / 4.0;
      z0 = (z0 + lam) / 4.0;
      mul *= 4.0;
    }
    final X = (A0 - x) / (mul * An);
    final Y = (A0 - y) / (mul * An);
    final Z = -(X + Y);
    final E2 = X * Y - Z * Z;
    final E3 = X * Y * Z;
    return (E3 * (6930.0 * E3 + E2 * (15015.0 * E2 - 16380.0) + 17160.0) +
            E2 * ((10010.0 - 5775.0 * E2) * E2 - 24024.0) +
            240240.0) /
        (240240.0 * math.sqrt(An));
  }

  /// RF(x, y) – degenerate form, Carlson eqs 2.36–2.38.
  static double _RF2(double x, double y) {
    final tolRG0 = 2.7 * math.sqrt(2.220446049250313e-16 * 0.01);
    double xn = math.sqrt(x), yn = math.sqrt(y);
    if (xn < yn) {
      final tmp = xn; xn = yn; yn = tmp;
    }
    while ((xn - yn).abs() > tolRG0 * xn) {
      final t = (xn + yn) / 2.0;
      yn = math.sqrt(xn * yn);
      xn = t;
    }
    return math.pi / (xn + yn);
  }

  /// RG(x, y) – two-argument version.
  static double _RG(double x, double y) {
    final tolRG0 = 2.7 * math.sqrt(2.220446049250313e-16 * 0.01);
    double x0 = math.sqrt(math.max(x, y));
    double y0 = math.sqrt(math.min(x, y));
    double xn = x0, yn = y0;
    double s = 0.0, mul = 0.25;
    while ((xn - yn).abs() > tolRG0 * xn) {
      final t = (xn + yn) / 2.0;
      yn = math.sqrt(xn * yn);
      xn = t;
      mul *= 2.0;
      final tDiff = xn - yn;
      s += mul * tDiff * tDiff;
    }
    return ((x0 + y0) / 2.0 * ((x0 + y0) / 2.0) - s) * math.pi / (2.0 * (xn + yn));
  }

  /// RD(x, y, z) – Carlson eqs 2.28–2.34.
  static double _RD(double x, double y, double z) {
    double A0 = (x + y + 3.0 * z) / 5.0;
    double An = A0;
    double Q = math.max(math.max((A0 - x).abs(), (A0 - y).abs()), (A0 - z).abs()) / _tolRD;
    double x0 = x, y0 = y, z0 = z;
    double mul = 1.0;
    double s = 0.0;
    while (Q >= mul * An.abs()) {
      final sq0 = math.sqrt(x0);
      final sq1 = math.sqrt(y0);
      final sq2 = math.sqrt(z0);
      final lam = sq0 * sq1 + sq1 * sq2 + sq2 * sq0;
      s += 1.0 / (mul * sq2 * (z0 + lam));
      An = (An + lam) / 4.0;
      x0 = (x0 + lam) / 4.0;
      y0 = (y0 + lam) / 4.0;
      z0 = (z0 + lam) / 4.0;
      mul *= 4.0;
    }
    final X = (A0 - x) / (mul * An);
    final Y = (A0 - y) / (mul * An);
    final Z = -(X + Y) / 3.0;
    final E2 = X * Y - 6.0 * Z * Z;
    final E3 = (3.0 * X * Y - 8.0 * Z * Z) * Z;
    final E4 = 3.0 * (X * Y - Z * Z) * Z * Z;
    final E5 = X * Y * Z * Z * Z;
    return ((471240.0 - 540540.0 * E2) * E5 +
            (612612.0 * E2 - 540540.0 * E3 - 556920.0) * E4 +
            E3 * (306306.0 * E3 + E2 * (675675.0 * E2 - 706860.0) + 680680.0) +
            E2 * ((417690.0 - 255255.0 * E2) * E2 - 875160.0) +
            4084080.0) /
        (4084080.0 * mul * An * math.sqrt(An)) +
        3.0 * s;
  }

  // ---------------------------------------------------------------------------
  // Small utility
  // ---------------------------------------------------------------------------
  static double _sinh(double x) {
    final e = math.exp(x);
    return (e - 1.0 / e) / 2.0;
  }

  static double _cosh(double x) {
    final e = math.exp(x);
    return (e + 1.0 / e) / 2.0;
  }

  static double _tanh(double x) {
    final e = math.exp(2.0 * x);
    return (e - 1.0) / (e + 1.0);
  }

  static double _copySign(double x, double y) {
    return x.abs() * (y.isNegative ? -1.0 : 1.0);
  }
}
