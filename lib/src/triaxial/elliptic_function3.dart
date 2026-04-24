// elliptic_function3.dart
// Dart port of GeographicLib::EllipticFunction from EllipticFunction.{hpp,cpp}
// and EllipticFunction.js (triaxial version)
//
// ==========================================================================
// Original C++/JS: Copyright (c) Charles Karney (2008-2024)
//   <karney@alum.mit.edu>, MIT/X11 License.
//   https://geographiclib.sourceforge.io/
//
// Dart Port:
//   Copyright (c) 2026 Friedhold Matz <fmatz.com@gmail.com>
//   Ported to Dart in 2026.
// Licence:
//   This Dart port is licensed under the MIT/X11 License (same as original).
//
// Algorithms:
//   Carlson (1995): doi:10.1007/BF02198293 — symmetric elliptic integrals
//   Bulirsch (1965): doi:10.1007/BF01397975 — sncndn (AGM method)
//   DLMF: https://dlmf.nist.gov/19, https://dlmf.nist.gov/22
// ==========================================================================

import 'dart:math' as math;

// ignore_for_file: non_constant_identifier_names

// ---------------------------------------------------------------------------
// Module-level constants
// ---------------------------------------------------------------------------

const double _pi = math.pi;
const double _machEps = 2.220446049250313e-16; // machine epsilon
const double _eps = _machEps; // legacy alias (used by top-level Carlson functions)
const double _inf = double.infinity;
const int _num = 25; // AGM depth (Bulirsch)

// ---------------------------------------------------------------------------
// Hyperbolic helpers (dart:math has no sinh/cosh/tanh)
// ---------------------------------------------------------------------------

double sinh3(double x) => (math.exp(x) - math.exp(-x)) / 2.0;
double cosh3(double x) => (math.exp(x) + math.exp(-x)) / 2.0;
double tanh3(double x) {
  if (x.isInfinite) return x.sign.toDouble();
  final e2 = math.exp(2.0 * x);
  return (e2 - 1.0) / (e2 + 1.0);
}

// ---------------------------------------------------------------------------
// Carlson symmetric integrals (static / top-level)
// ---------------------------------------------------------------------------

/// RF(x, y, z) — Carlson symmetric integral of the first kind.
double rf(double x, double y, double z) {
  final tolRF = math.pow(3 * _machEps * 0.01, 1.0 / 8.0) as double;
  final A0 = (x + y + z) / 3.0;
  var An = A0;
  final Q = math.max(
          math.max((A0 - x).abs(), (A0 - y).abs()), (A0 - z).abs()) /
      tolRF;
  var x0 = x, y0 = y, z0 = z, mul = 1.0;
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
  return (E3 * (6930 * E3 + E2 * (15015 * E2 - 16380) + 17160) +
          E2 * ((10010 - 5775 * E2) * E2 - 24024) +
          240240) /
      (240240 * math.sqrt(An));
}

/// RF(x, y) — degenerate case z=0 (Carlson eq. 2.36–2.38).
double rf2(double x, double y) {
  final tolRG0 = 2.7 * math.sqrt(_machEps * 0.01);
  var xn = math.sqrt(x), yn = math.sqrt(y);
  if (xn < yn) {
    final tmp = xn;
    xn = yn;
    yn = tmp;
  }
  while ((xn - yn).abs() > tolRG0 * xn) {
    final t = (xn + yn) / 2.0;
    yn = math.sqrt(xn * yn);
    xn = t;
  }
  return _pi / (xn + yn);
}

/// RC(x, y).
double rc(double x, double y) {
  if (!(x >= y)) {
    return math.atan(math.sqrt((y - x) / x)) / math.sqrt(y - x);
  }
  if (x == y) return 1.0 / math.sqrt(y);
  final val = y > 0
      ? math.sqrt((x - y) / y)
      : math.sqrt(-x / y);
  // asinh(val) = log(val + sqrt(val*val + 1))
  return math.log(val + math.sqrt(val * val + 1.0)) / math.sqrt(x - y);
}

/// RG(x, y, z) — Carlson symmetric integral of the second kind.
double rg(double x, double y, double z) {
  if (x == 0.0) return rg2(y, z);
  if (y == 0.0) return rg2(z, x);
  if (z == 0.0) return rg2(x, y);
  return (z * rf(x, y, z) - (x - z) * (y - z) * rd(x, y, z) / 3.0 +
          math.sqrt(x * y / z)) /
      2.0;
}

/// RG(x, y) — degenerate case (Carlson eq. 2.36–2.39).
double rg2(double x, double y) {
  final tolRG0 = 2.7 * math.sqrt(_machEps * 0.01);
  var xn = math.sqrt(math.max(x, y));
  var yn = math.sqrt(math.min(x, y));
  var s = 0.0, mul = 0.25;
  while ((xn - yn).abs() > tolRG0 * xn) {
    final t = (xn + yn) / 2.0;
    yn = math.sqrt(xn * yn);
    xn = t;
    mul *= 2.0;
    final diff = xn - yn;
    s += mul * diff * diff;
  }
  final avg = (math.sqrt(x) + math.sqrt(y)) / 2.0;
  return (avg * avg - s) * _pi / (2.0 * (xn + yn));
}

/// RJ(x, y, z, p) — Carlson symmetric integral of the third kind.
double rj(double x, double y, double z, double p) {
  final tolRD = math.pow(0.2 * _machEps * 0.01, 1.0 / 8.0) as double;
  final A0 = (x + y + z + 2.0 * p) / 5.0;
  var An = A0;
  final delta = (p - x) * (p - y) * (p - z);
  final Q = math.max(
          math.max(
              math.max((A0 - x).abs(), (A0 - y).abs()), (A0 - z).abs()),
          (A0 - p).abs()) /
      tolRD;
  var x0 = x, y0 = y, z0 = z, p0 = p;
  var mul = 1.0, mul3 = 1.0, s = 0.0;
  while (Q >= mul * An.abs()) {
    final lam = math.sqrt(x0) * math.sqrt(y0) +
        math.sqrt(y0) * math.sqrt(z0) +
        math.sqrt(z0) * math.sqrt(x0);
    final d0 = (math.sqrt(p0) + math.sqrt(x0)) *
        (math.sqrt(p0) + math.sqrt(y0)) *
        (math.sqrt(p0) + math.sqrt(z0));
    final e0 = delta / (mul3 * d0 * d0);
    s += rc(1.0, 1.0 + e0) / (mul * d0);
    An = (An + lam) / 4.0;
    x0 = (x0 + lam) / 4.0;
    y0 = (y0 + lam) / 4.0;
    z0 = (z0 + lam) / 4.0;
    p0 = (p0 + lam) / 4.0;
    mul *= 4.0;
    mul3 *= 64.0;
  }
  final X = (A0 - x) / (mul * An);
  final Y = (A0 - y) / (mul * An);
  final Z = (A0 - z) / (mul * An);
  final P = -(X + Y + Z) / 2.0;
  final E2 = X * Y + X * Z + Y * Z - 3.0 * P * P;
  final E3 = X * Y * Z + 2.0 * P * (E2 + 2.0 * P * P);
  final E4 = (2.0 * X * Y * Z + P * (E2 + 3.0 * P * P)) * P;
  final E5 = X * Y * Z * P * P;
  return ((471240 - 540540 * E2) * E5 +
          (612612 * E2 - 540540 * E3 - 556920) * E4 +
          E3 * (306306 * E3 + E2 * (675675 * E2 - 706860) + 680680) +
          E2 * ((417690 - 255255 * E2) * E2 - 875160) +
          4084080) /
      (4084080 * mul * An * math.sqrt(An)) +
      6.0 * s;
}

/// RD(x, y, z) = RJ(x, y, z, z) (symmetric integral of the second kind).
double rd(double x, double y, double z) {
  final tolRD = math.pow(0.2 * _eps * 0.01, 1.0 / 8.0) as double;
  final A0 = (x + y + 3.0 * z) / 5.0;
  var An = A0;
  final Q = math.max(
          math.max((A0 - x).abs(), (A0 - y).abs()), (A0 - z).abs()) /
      tolRD;
  var x0 = x, y0 = y, z0 = z, mul = 1.0, s = 0.0;
  while (Q >= mul * An.abs()) {
    final lam = math.sqrt(x0) * math.sqrt(y0) +
        math.sqrt(y0) * math.sqrt(z0) +
        math.sqrt(z0) * math.sqrt(x0);
    s += 1.0 / (mul * math.sqrt(z0) * (z0 + lam));
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
  return ((471240 - 540540 * E2) * E5 +
          (612612 * E2 - 540540 * E3 - 556920) * E4 +
          E3 * (306306 * E3 + E2 * (675675 * E2 - 706860) + 680680) +
          E2 * ((417690 - 255255 * E2) * E2 - 875160) +
          4084080) /
      (4084080 * mul * An * math.sqrt(An)) +
      3.0 * s;
}

// ---------------------------------------------------------------------------
// EllipticFunction class (triaxial version)
// ---------------------------------------------------------------------------

/// Elliptic integrals and Jacobi elliptic functions for the triaxial geodesic.
class EllipticFunction3 {
  double _k2;
  double _kp2;
  double _alpha2;
  double _alphap2;
  double _eps;

  // Complete integrals
  double _kKc = 0.0;
  double _eEc = 0.0;
  double _dDc = 0.0;
  double _pPic = 0.0;
  double _gGc = 0.0;
  double _hHc = 0.0;

  /// Construct with [k2] and optional [alpha2].
  EllipticFunction3([double k2 = 0.0, double alpha2 = 0.0])
      : _k2 = k2,
        _kp2 = 1.0 - k2,
        _alpha2 = alpha2,
        _alphap2 = 1.0 - alpha2,
        _eps = 0.0 {
    _reset4(k2, alpha2, 1.0 - k2, 1.0 - alpha2);
  }

  /// Full precision reset (avoids cancellation for small kp2).
  EllipticFunction3.full(double k2, double alpha2, double kp2, double alphap2)
      : _k2 = k2,
        _kp2 = kp2,
        _alpha2 = alpha2,
        _alphap2 = alphap2,
        _eps = 0.0 {
    _reset4(k2, alpha2, kp2, alphap2);
  }

  // -------------------------------------------------------------------------
  // Reset
  // -------------------------------------------------------------------------

  /// Reset with two arguments.
  void reset([double k2 = 0.0, double alpha2 = 0.0]) {
    _reset4(k2, alpha2, 1.0 - k2, 1.0 - alpha2);
  }

  /// Reset with four arguments (full precision for small kp2).
  void reset4(double k2, double alpha2, double kp2, double alphap2) {
    _reset4(k2, alpha2, kp2, alphap2);
  }

  void _reset4(double k2, double alpha2, double kp2, double alphap2) {
    _k2 = k2;
    _kp2 = kp2;
    _alpha2 = alpha2;
    _alphap2 = alphap2;
    _eps = k2 / (math.sqrt(kp2) + 1.0) * (math.sqrt(kp2) + 1.0);

    if (k2 != 0.0) {
      _kKc = kp2 != 0.0 ? rf2(kp2, 1.0) : _inf;
      _eEc = kp2 != 0.0 ? 2.0 * rg2(kp2, 1.0) : 1.0;
      _dDc = kp2 != 0.0 ? rd(0.0, kp2, 1.0) / 3.0 : _inf;
    } else {
      _kKc = _eEc = _pi / 2.0;
      _dDc = _pi / 4.0;
    }

    if (alpha2 != 0.0) {
      final rjVal = (kp2 != 0.0 && alphap2 != 0.0)
          ? rj(0.0, kp2, 1.0, alphap2)
          : _inf;
      final rcVal = kp2 != 0.0
          ? 0.0
          : (alphap2 != 0.0 ? rc(1.0, alphap2) : _inf);
      _pPic = kp2 != 0.0 ? _kKc + alpha2 * rjVal / 3.0 : _inf;
      _gGc = kp2 != 0.0 ? _kKc + (alpha2 - k2) * rjVal / 3.0 : rcVal;
      _hHc = kp2 != 0.0
          ? _kKc - (alphap2 != 0.0 ? alphap2 * rjVal : 0.0) / 3.0
          : rcVal;
    } else {
      _pPic = _kKc;
      _gGc = _eEc;
      _hHc = kp2 == 1.0
          ? _pi / 4.0
          : (kp2 == 0.0 ? 1.0 : kp2 * rd(0.0, 1.0, kp2) / 3.0);
    }
  }

  // -------------------------------------------------------------------------
  // Inspectors
  // -------------------------------------------------------------------------

  double get k2 => _k2;
  double get kp2 => _kp2;
  double get alpha2 => _alpha2;
  double get alphap2 => _alphap2;

  // Complete integrals
  double get K => _kKc;
  double get E => _eEc;
  double get D => _dDc;
  double get KE => _k2 * _dDc;
  double get Pi => _pPic;
  double get G => _gGc;
  double get H => _hHc;

  // -------------------------------------------------------------------------
  // Delta amplitude
  // -------------------------------------------------------------------------

  /// Δ(sn, cn) = sqrt(1 - k²·sin²φ).
  double delta(double sn, double cn) {
    return math.sqrt(_k2 < 0
        ? 1.0 - _k2 * sn * sn
        : _kp2 + _k2 * cn * cn);
  }

  // -------------------------------------------------------------------------
  // Jacobi elliptic functions: sncndn  (Bulirsch AGM method)
  // -------------------------------------------------------------------------

  /// Compute sn(x,k), cn(x,k), dn(x,k).
  ({double sn, double cn, double dn}) sncndn(double x) {
    final tolJAC = math.sqrt(_machEps * 0.01);
    var kp2 = _kp2;
    var d = 0.0;
    if (kp2 != 0.0) {
      var mc = kp2;
      if (mc < 0) {
        d = 1.0 - mc;
        mc /= -d;
        d = math.sqrt(d);
        x *= d;
      }
      final m = List<double>.filled(_num, 0.0);
      final n = List<double>.filled(_num, 0.0);
      var l = 0;
      var a = 1.0, c = 0.0;
      for (; l < _num; ++l) {
        m[l] = a;
        n[l] = mc = math.sqrt(mc);
        c = (a + mc) / 2.0;
        if ((a - mc).abs() <= tolJAC * a) {
          ++l;
          break;
        }
        mc *= a;
        a = c;
      }
      x *= c;
      double sn = math.sin(x), cn = math.cos(x), dn = 1.0;
      if (sn != 0.0) {
        a = cn / sn;
        c *= a;
        while (l-- > 0) {
          final b = m[l];
          a *= c;
          c *= dn;
          dn = (n[l] + a) / (b + a);
          a = c / b;
        }
        a = 1.0 / math.sqrt(c * c + 1.0);
        sn = sn < 0 ? -a : a;
        cn = c * sn;
        if (kp2 < 0) {
          final tmp = cn;
          cn = dn;
          dn = tmp;
          sn /= d;
        }
      }
      return (sn: sn, cn: cn, dn: dn);
    } else {
      final t = 1.0 / cosh3(x);
      return (sn: tanh3(x), cn: t, dn: t);
    }
  }

  // -------------------------------------------------------------------------
  // Jacobi amplitude am(x) — DLMF 22.20(ii)
  // -------------------------------------------------------------------------

  /// am(x, k): Jacobi amplitude function.
  double am(double x) {
    final tolJAC = math.pow(_machEps, 0.75) as double;
    var k2 = _k2, kp2 = _kp2;
    if (k2 == 0.0) return x;
    if (kp2 == 0.0) return math.atan(sinh3(x));
    if (k2 < 0) {
      k2 = -k2 / kp2;
      kp2 = 1.0 / kp2;
      x *= math.sqrt(_kp2);
    }
    final a = List<double>.filled(_num + 1, 0.0);
    var b = math.sqrt(kp2);
    final c = List<double>.filled(_num + 1, 0.0);
    a[0] = 1.0;
    c[0] = math.sqrt(k2);
    var l = 1;
    for (; l < _num; ++l) {
      a[l] = (a[l - 1] + b) / 2.0;
      c[l] = (a[l - 1] - b) / 2.0;
      b = math.sqrt(a[l - 1] * b);
      if (!(c[l] > tolJAC * a[l])) break;
    }
    var phi = a[l] * x * math.pow(2, l);
    var phi1 = 0.0;
    for (; l > 0; --l) {
      phi1 = phi;
      phi = (phi + math.asin(c[l] * math.sin(phi) / a[l])) / 2.0;
    }
    return _k2 < 0 ? phi1 - phi : phi;
  }

  /// am(x) with sn, cn, dn output.
  ({double phi, double sn, double cn, double dn}) amAll(double x) {
    final phi = am(x);
    double sn, cn, dn;
    if (_kp2 == 0.0) {
      sn = tanh3(x);
      cn = dn = 1.0 / cosh3(x);
    } else {
      sn = math.sin(phi);
      cn = math.cos(phi);
      dn = delta(sn, cn);
    }
    return (phi: phi, sn: sn, cn: cn, dn: dn);
  }

  // -------------------------------------------------------------------------
  // Incomplete integrals in Jacobi form
  // -------------------------------------------------------------------------

  double fSncd(double sn, double cn, double dn) {
    final cn2 = cn * cn, dn2 = dn * dn;
    var fi = cn2 != 0.0 ? sn.abs() * rf(cn2, dn2, 1.0) : _kKc;
    if (cn < 0) fi = 2.0 * _kKc - fi;
    return sn < 0 ? -fi : fi;
  }

  double eSncd(double sn, double cn, double dn) {
    final k2 = _k2, kp2 = _kp2;
    final cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn;
    double ei;
    if (cn2 != 0.0) {
      final asn = sn.abs();
      if (k2 <= 0) {
        ei = asn * (rf(cn2, dn2, 1.0) - k2 * sn2 * rd(cn2, dn2, 1.0) / 3.0);
      } else if (kp2 >= 0) {
        ei = asn * (kp2 * rf(cn2, dn2, 1.0) +
            k2 * kp2 * sn2 * rd(cn2, 1.0, dn2) / 3.0 +
            k2 * cn.abs() / dn);
      } else {
        ei = asn * (-kp2 * sn2 * rd(dn2, 1.0, cn2) / 3.0 + dn / cn.abs());
      }
    } else {
      ei = _eEc;
    }
    if (cn < 0) ei = 2.0 * _eEc - ei;
    return sn < 0 ? -ei : ei;
  }

  double dSncd(double sn, double cn, double dn) {
    final cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn;
    var di = cn2 != 0.0
        ? sn.abs() * sn2 * rd(cn2, dn2, 1.0) / 3.0
        : _dDc;
    if (cn < 0) di = 2.0 * _dDc - di;
    return sn < 0 ? -di : di;
  }

  double piSncd(double sn, double cn, double dn) {
    final alpha2 = _alpha2, alphap2 = _alphap2;
    final cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn;
    var pii = cn2 != 0.0
        ? sn.abs() *
            (rf(cn2, dn2, 1.0) +
                alpha2 * sn2 * rj(cn2, dn2, 1.0, cn2 + alphap2 * sn2) / 3.0)
        : _pPic;
    if (cn < 0) pii = 2.0 * _pPic - pii;
    return sn < 0 ? -pii : pii;
  }

  double gSncd(double sn, double cn, double dn) {
    final k2 = _k2, alpha2 = _alpha2, alphap2 = _alphap2;
    final cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn;
    var gi = cn2 != 0.0
        ? sn.abs() *
            (rf(cn2, dn2, 1.0) +
                (alpha2 - k2) *
                    sn2 *
                    rj(cn2, dn2, 1.0, cn2 + alphap2 * sn2) /
                    3.0)
        : _gGc;
    if (cn < 0) gi = 2.0 * _gGc - gi;
    return sn < 0 ? -gi : gi;
  }

  double hSncd(double sn, double cn, double dn) {
    final alphap2 = _alphap2;
    final cn2 = cn * cn, dn2 = dn * dn, sn2 = sn * sn;
    var hi = cn2 != 0.0
        ? sn.abs() *
            (rf(cn2, dn2, 1.0) -
                alphap2 * sn2 * rj(cn2, dn2, 1.0, cn2 + alphap2 * sn2) / 3.0)
        : _hHc;
    if (cn < 0) hi = 2.0 * _hHc - hi;
    return sn < 0 ? -hi : hi;
  }

  // -------------------------------------------------------------------------
  // Periodic (delta) variants
  // -------------------------------------------------------------------------

  double deltaF(double sn, double cn, double dn) {
    if (cn < 0) {
      cn = -cn;
      sn = -sn;
    }
    return fSncd(sn, cn, dn) * (_pi / 2.0) / _kKc - math.atan2(sn, cn);
  }

  double deltaE(double sn, double cn, double dn) {
    if (cn < 0) {
      cn = -cn;
      sn = -sn;
    }
    return eSncd(sn, cn, dn) * (_pi / 2.0) / _eEc - math.atan2(sn, cn);
  }

  double deltaPi(double sn, double cn, double dn) {
    if (cn < 0) {
      cn = -cn;
      sn = -sn;
    }
    return piSncd(sn, cn, dn) * (_pi / 2.0) / _pPic - math.atan2(sn, cn);
  }

  double deltaD(double sn, double cn, double dn) {
    if (cn < 0) {
      cn = -cn;
      sn = -sn;
    }
    return dSncd(sn, cn, dn) * (_pi / 2.0) / _dDc - math.atan2(sn, cn);
  }

  double deltaG(double sn, double cn, double dn) {
    if (cn < 0) {
      cn = -cn;
      sn = -sn;
    }
    return gSncd(sn, cn, dn) * (_pi / 2.0) / _gGc - math.atan2(sn, cn);
  }

  double deltaH(double sn, double cn, double dn) {
    if (cn < 0) {
      cn = -cn;
      sn = -sn;
    }
    return hSncd(sn, cn, dn) * (_pi / 2.0) / _hHc - math.atan2(sn, cn);
  }

  // -------------------------------------------------------------------------
  // Incomplete integrals in angle form
  // -------------------------------------------------------------------------

  double fPhi(double phi) {
    if (_k2 == 0.0) return phi;
    if (_kp2 == 0.0) return math.log(math.tan(phi) + 1.0 / math.cos(phi));
    final sn = math.sin(phi), cn = math.cos(phi), dn = delta(sn, cn);
    return phi.abs() < _pi
        ? fSncd(sn, cn, dn)
        : (deltaF(sn, cn, dn) + phi) * _kKc / (_pi / 2.0);
  }

  double ePhi(double phi) {
    if (_k2 == 0.0) return phi;
    final sn = math.sin(phi), cn = math.cos(phi), dn = delta(sn, cn);
    return phi.abs() < _pi
        ? eSncd(sn, cn, dn)
        : (deltaE(sn, cn, dn) + phi) * _eEc / (_pi / 2.0);
  }

  double edPhi(double ang) {
    final n = ((ang - ((ang % 360.0) +
            (ang % 360.0 < -180.0
                ? 360.0
                : (ang % 360.0 > 180.0 ? -360.0 : 0.0)))) /
            360.0)
        .roundToDouble();
    final sn = math.sin(ang * _pi / 180.0);
    final cn = math.cos(ang * _pi / 180.0);
    return eSncd(sn, cn, delta(sn, cn)) + 4.0 * _eEc * n;
  }

  double piPhi(double phi) {
    final sn = math.sin(phi), cn = math.cos(phi), dn = delta(sn, cn);
    return phi.abs() < _pi
        ? piSncd(sn, cn, dn)
        : (deltaPi(sn, cn, dn) + phi) * _pPic / (_pi / 2.0);
  }

  double dPhi(double phi) {
    final sn = math.sin(phi), cn = math.cos(phi), dn = delta(sn, cn);
    return phi.abs() < _pi
        ? dSncd(sn, cn, dn)
        : (deltaD(sn, cn, dn) + phi) * _dDc / (_pi / 2.0);
  }

  double gPhi(double phi) {
    final sn = math.sin(phi), cn = math.cos(phi), dn = delta(sn, cn);
    return phi.abs() < _pi
        ? gSncd(sn, cn, dn)
        : (deltaG(sn, cn, dn) + phi) * _gGc / (_pi / 2.0);
  }

  double hPhi(double phi) {
    final sn = math.sin(phi), cn = math.cos(phi), dn = delta(sn, cn);
    return phi.abs() < _pi
        ? hSncd(sn, cn, dn)
        : (deltaH(sn, cn, dn) + phi) * _hHc / (_pi / 2.0);
  }

  // -------------------------------------------------------------------------
  // Einv: inverse of incomplete E
  // -------------------------------------------------------------------------

  double einv(double x) {
    final tolJAC = math.sqrt(_machEps / 100.0);
    final ec = _eEc;
    final n = (x / (2.0 * ec) + 0.5).floorToDouble();
    x -= 2.0 * ec * n;
    var phi = _pi * x / (2.0 * ec);
    phi -= _eps * math.sin(2.0 * phi) / 2.0;
    for (var i = 0; i < _num; ++i) {
      final sn = math.sin(phi), cn = math.cos(phi), dn = delta(sn, cn);
      final err = (eSncd(sn, cn, dn) - x) / dn;
      phi -= err;
      if (err.abs() <= tolJAC) break;
    }
    return n * _pi + phi;
  }

  double deltaEinv(double stau, double ctau) {
    if (ctau < 0) {
      ctau = -ctau;
      stau = -stau;
    }
    final tau = math.atan2(stau, ctau);
    return einv(tau * _eEc / (_pi / 2.0)) - tau;
  }
}




