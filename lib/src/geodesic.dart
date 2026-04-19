// geodesic.dart
// Dart port of Geodesic.js + GeodesicLine.js + PolygonArea.js
//   from geographiclib-geodesic v2.2.0
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
// ==========================================================================

import 'dart:math' as math;
import 'geo_math.dart';

// ignore_for_file: non_constant_identifier_names

const int _order = 6;
const int _nA1 = _order;
const int _nA2 = _order;
const int _nA3 = _order;
const int _nA3x = _order;
const int _nC3x = (_order * (_order - 1)) ~/ 2;
const int _nC4x = (_order * (_order + 1)) ~/ 2;
const int _maxit1 = 20;
const int _maxit2 = _maxit1 + GeoMath.digits + 10;
const double _tol0 = GeoMath.epsilon;
const double _tol1 = 200.0 * _tol0;
final double _tol2 = math.sqrt(_tol0);
const double _tolb = _tol0;
final double _xthresh = 1000.0 * math.sqrt(_tol0);

// Capability flags
const int capC1   = 1 << 0;
const int capC1p  = 1 << 1;
const int capC2   = 1 << 2;
const int capC3   = 1 << 3;
const int capC4   = 1 << 4;
const int capNone = 0;
const int capAll  = 0x1F;

// Output flags
const int gNone          = 0;
const int gArc           = 1 << 6;
const int gLatitude      = 1 << 7;
const int gLongitude     = (1 << 8) | capC3;
const int gAzimuth       = 1 << 9;
const int gDistance      = (1 << 10) | capC1;
const int gStandard      = gLatitude | gLongitude | gAzimuth | gDistance;
const int gDistanceIn    = (1 << 11) | capC1 | capC1p;
const int gReducedLength = (1 << 12) | capC1 | capC2;
const int gGeodesicScale = (1 << 13) | capC1 | capC2;
const int gArea          = (1 << 14) | capC4;
const int gAll           = 0x7F80 | capAll;
const int gLongUnroll    = 1 << 15;
const int gOutMask       = 0x7F80 | gLongUnroll;

// tiny_ = sqrt(MIN_POSITIVE / epsilon)
final double _tiny = math.sqrt(5e-324 / GeoMath.epsilon);

const int _nC1 = _order;
const int _nC1p = _order;
const int _nC2 = _order;
const int _nC3 = _order;
const int _nC4 = _order;

// ── Coefficient arrays ────────────────────────────────────────────────────────

final List<double> _A1m1fCoeff = [1, 4, 64, 0, 256].map((e) => e.toDouble()).toList();

double _A1m1f(double eps) {
  const p = _nA1 ~/ 2;
  final t = GeoMath.polyval(p, _A1m1fCoeff, 0, GeoMath.sq(eps)) / _A1m1fCoeff[p + 1];
  return (t + eps) / (1.0 - eps);
}

final List<double> _C1fCoeff = [
  -1, 6, -16, 32,
  -9, 64, -128, 2048,
  9, -16, 768,
  3, -5, 512,
  -7, 1280,
  -7, 2048,
].map((e) => e.toDouble()).toList();

void _C1f(double eps, List<double> c) {
  final eps2 = GeoMath.sq(eps);
  double d = eps;
  int o = 0;
  for (int l = 1; l <= _nC1; ++l) {
    final p = (_nC1 - l) ~/ 2;
    c[l] = d * GeoMath.polyval(p, _C1fCoeff, o, eps2) / _C1fCoeff[o + p + 1];
    o += p + 2;
    d *= eps;
  }
}

final List<double> _C1pfCoeff = [
  205, -432, 768, 1536,
  4005, -4736, 3840, 12288,
  -225, 116, 384,
  -7173, 2695, 7680,
  3467, 7680,
  38081, 61440,
].map((e) => e.toDouble()).toList();

void _C1pf(double eps, List<double> c) {
  final eps2 = GeoMath.sq(eps);
  double d = eps;
  int o = 0;
  for (int l = 1; l <= _nC1p; ++l) {
    final p = (_nC1p - l) ~/ 2;
    c[l] = d * GeoMath.polyval(p, _C1pfCoeff, o, eps2) / _C1pfCoeff[o + p + 1];
    o += p + 2;
    d *= eps;
  }
}

final List<double> _A2m1fCoeff = [-11, -28, -192, 0, 256].map((e) => e.toDouble()).toList();

double _A2m1f(double eps) {
  const p = _nA2 ~/ 2;
  final t = GeoMath.polyval(p, _A2m1fCoeff, 0, GeoMath.sq(eps)) / _A2m1fCoeff[p + 1];
  return (t - eps) / (1.0 + eps);
}

final List<double> _C2fCoeff = [
  1, 2, 16, 32,
  35, 64, 384, 2048,
  15, 80, 768,
  7, 35, 512,
  63, 1280,
  77, 2048,
].map((e) => e.toDouble()).toList();

void _C2f(double eps, List<double> c) {
  final eps2 = GeoMath.sq(eps);
  double d = eps;
  int o = 0;
  for (int l = 1; l <= _nC2; ++l) {
    final p = (_nC2 - l) ~/ 2;
    c[l] = d * GeoMath.polyval(p, _C2fCoeff, o, eps2) / _C2fCoeff[o + p + 1];
    o += p + 2;
    d *= eps;
  }
}

// SinCosSeries: Clenshaw summation
double sinCosSeries(bool sinp, double sinx, double cosx, List<double> c) {
  int k = c.length;
  int n = k - (sinp ? 1 : 0);
  final ar = 2.0 * (cosx - sinx) * (cosx + sinx); // 2*cos(2x)
  double y0 = (n & 1) != 0 ? c[--k] : 0.0;
  double y1 = 0.0;
  n = n ~/ 2;
  while (n-- > 0) {
    y1 = ar * y0 - y1 + c[--k];
    y0 = ar * y1 - y0 + c[--k];
  }
  return sinp ? 2.0 * sinx * cosx * y0 : cosx * (y0 - y1);
}

double _astroid(double x, double y) {
  final p = GeoMath.sq(x);
  final q = GeoMath.sq(y);
  final r = (p + q - 1.0) / 6.0;
  if (!(q == 0.0 && r <= 0.0)) {
    final S = p * q / 4.0;
    final r2 = GeoMath.sq(r);
    final r3 = r * r2;
    final disc = S * (S + 2.0 * r3);
    double u = r;
    if (disc >= 0.0) {
      double T3 = S + r3;
      T3 += T3 < 0.0 ? -math.sqrt(disc) : math.sqrt(disc);
      final T = GeoMath.cbrt(T3);
      u += T + (T != 0.0 ? r2 / T : 0.0);
    } else {
      final ang = math.atan2(math.sqrt(-disc), -(S + r3));
      u += 2.0 * r * math.cos(ang / 3.0);
    }
    final v = math.sqrt(GeoMath.sq(u) + q);
    final uv = u < 0.0 ? q / (v - u) : u + v;
    final w = (uv - q) / (2.0 * v);
    return uv / (math.sqrt(uv + GeoMath.sq(w)) + w);
  } else {
    return 0.0;
  }
}

/// Result of a geodesic inverse or direct computation.
class GeodesicData {
  double lat1, lon1, azi1;
  double lat2, lon2, azi2;
  double s12, a12;
  double m12, M12, M21, S12;

  GeodesicData()
      : lat1 = double.nan,
        lon1 = double.nan,
        azi1 = double.nan,
        lat2 = double.nan,
        lon2 = double.nan,
        azi2 = double.nan,
        s12 = double.nan,
        a12 = double.nan,
        m12 = double.nan,
        M12 = double.nan,
        M21 = double.nan,
        S12 = double.nan;
}

// Internal struct used by InverseInt
class _InverseResult {
  GeodesicData vals;
  double salp1, calp1, salp2, calp2;
  _InverseResult(this.vals, this.salp1, this.calp1, this.salp2, this.calp2);
}

/// The A3 coefficient array; populated in Geodesic constructor.
final List<double> _A3coeff_global = [
  -3, 128,
  -2, -3, 64,
  -1, -3, -1, 16,
  3, -1, -2, 8,
  1, -1, 2,
  1, 1,
].map((e) => e.toDouble()).toList();

final List<double> _C3coeff_global = [
  3, 128,
  2, 5, 128,
  -1, 3, 3, 64,
  -1, 0, 1, 8,
  -1, 1, 4,
  5, 256,
  1, 3, 128,
  -3, -2, 3, 64,
  1, -3, 2, 32,
  7, 512,
  -10, 9, 384,
  5, -9, 5, 192,
  7, 512,
  -14, 7, 512,
  21, 2560,
].map((e) => e.toDouble()).toList();

final List<double> _C4coeff_global = [
  97, 15015,
  1088, 156, 45045,
  -224, -4784, 1573, 45045,
  -10656, 14144, -4576, -858, 45045,
  64, 624, -4576, 6864, -3003, 15015,
  100, 208, 572, 3432, -12012, 30030, 45045,
  1, 9009,
  -2944, 468, 135135,
  5792, 1040, -1287, 135135,
  5952, -11648, 9152, -2574, 135135,
  -64, -624, 4576, -6864, 3003, 135135,
  8, 10725,
  1856, -936, 225225,
  -8448, 4992, -1144, 225225,
  -1440, 4160, -4576, 1716, 225225,
  -136, 63063,
  1024, -208, 105105,
  3584, -3328, 1144, 315315,
  -128, 135135,
  -2560, 832, 405405,
  128, 99099,
].map((e) => e.toDouble()).toList();

/// Geodesic solver for an ellipsoid of revolution.
class Geodesic {
  /// Equatorial radius (meters).
  final double a;

  /// Flattening.
  final double f;

  final double _f1;
  final double _e2;
  final double _ep2;
  final double _n;
  final double _b;
  final double _c2;
  final double _etol2;

  late final List<double> _A3x;
  late final List<double> _C3x;
  late final List<double> _C4x;

  /// Standard WGS84 instance.
  static final Geodesic WGS84 = Geodesic(6378137.0, 1.0 / 298.257223563);

  Geodesic(this.a, this.f)
      : _f1 = 1.0 - f,
        _e2 = f * (2.0 - f),
        _ep2 = f * (2.0 - f) / ((1.0 - f) * (1.0 - f)),
        _n = f / (2.0 - f),
        _b = a * (1.0 - f),
        _c2 = _computeC2(a, f),
        _etol2 = _computeEtol2(f) {
    if (!a.isFinite || a <= 0) throw ArgumentError('Equatorial radius is not positive');
    if (!_b.isFinite || _b <= 0) throw ArgumentError('Polar semi-axis is not positive');
    _A3x = List<double>.filled(_nA3x, 0.0);
    _C3x = List<double>.filled(_nC3x, 0.0);
    _C4x = List<double>.filled(_nC4x, 0.0);
    _computeA3coeff();
    _computeC3coeff();
    _computeC4coeff();
  }

  static double _computeC2(double a, double f) {
    final f1 = 1.0 - f;
    final b = a * f1;
    final e2 = f * (2.0 - f);
    double e2term;
    if (e2 == 0.0) {
      e2term = 1.0;
    } else if (e2 > 0.0) {
      e2term = GeoMath.atanh(math.sqrt(e2)) / math.sqrt(e2.abs());
    } else {
      e2term = math.atan(math.sqrt(-e2)) / math.sqrt(e2.abs());
    }
    return (GeoMath.sq(a) + GeoMath.sq(b) * e2term) / 2.0;
  }

  static double _computeEtol2(double f) {
    return 0.1 * math.sqrt(_tol0) /
        math.sqrt(math.max(0.001, f.abs()) * math.min(1.0, 1.0 - f / 2.0) / 2.0);
  }

  void _computeA3coeff() {
    final coeff = _A3coeff_global;
    int o = 0, k = 0;
    for (int j = _nA3 - 1; j >= 0; --j) {
      final p = math.min(_nA3 - j - 1, j);
      _A3x[k++] = GeoMath.polyval(p, coeff, o, _n) / coeff[o + p + 1];
      o += p + 2;
    }
  }

  void _computeC3coeff() {
    final coeff = _C3coeff_global;
    int o = 0, k = 0;
    for (int l = 1; l < _nC3; ++l) {
      for (int j = _nC3 - 1; j >= l; --j) {
        final p = math.min(_nC3 - j - 1, j);
        _C3x[k++] = GeoMath.polyval(p, coeff, o, _n) / coeff[o + p + 1];
        o += p + 2;
      }
    }
  }

  void _computeC4coeff() {
    final coeff = _C4coeff_global;
    int o = 0, k = 0;
    for (int l = 0; l < _nC4; ++l) {
      for (int j = _nC4 - 1; j >= l; --j) {
        final p = _nC4 - j - 1;
        _C4x[k++] = GeoMath.polyval(p, coeff, o, _n) / coeff[o + p + 1];
        o += p + 2;
      }
    }
  }

  double _A3f(double eps) => GeoMath.polyval(_nA3x - 1, _A3x, 0, eps);

  void _C3f(double eps, List<double> c) {
    double mult = 1.0;
    int o = 0;
    for (int l = 1; l < _nC3; ++l) {
      final p = _nC3 - l - 1;
      mult *= eps;
      c[l] = mult * GeoMath.polyval(p, _C3x, o, eps);
      o += p + 1;
    }
  }

  void _C4f(double eps, List<double> c) {
    double mult = 1.0;
    int o = 0;
    for (int l = 0; l < _nC4; ++l) {
      final p = _nC4 - l - 1;
      c[l] = mult * GeoMath.polyval(p, _C4x, o, eps);
      o += p + 1;
      mult *= eps;
    }
  }

  // ── Lengths ────────────────────────────────────────────────────────────────

  ({double s12b, double m12b, double m0, double M12, double M21}) _lengths(
      double eps, double sig12,
      double ssig1, double csig1, double dn1,
      double ssig2, double csig2, double dn2,
      double cbet1, double cbet2, int outmask,
      List<double> C1a, List<double> C2a) {
    outmask &= gOutMask & 0x7F80; // OUT_MASK without LONG_UNROLL
    double m0x = 0.0, J12 = 0.0, A1 = 0.0, A2 = 0.0;
    double s12b = double.nan, m12b = double.nan, m0 = double.nan;
    double M12 = double.nan, M21 = double.nan;

    if ((outmask & (gDistance | gReducedLength | gGeodesicScale)) != 0) {
      A1 = _A1m1f(eps);
      _C1f(eps, C1a);
      if ((outmask & (gReducedLength | gGeodesicScale)) != 0) {
        A2 = _A2m1f(eps);
        _C2f(eps, C2a);
        m0x = A1 - A2;
        A2 = 1.0 + A2;
      }
      A1 = 1.0 + A1;
    }

    if ((outmask & gDistance) != 0) {
      final B1 = sinCosSeries(true, ssig2, csig2, C1a) -
          sinCosSeries(true, ssig1, csig1, C1a);
      s12b = A1 * (sig12 + B1);
      if ((outmask & (gReducedLength | gGeodesicScale)) != 0) {
        final B2 = sinCosSeries(true, ssig2, csig2, C2a) -
            sinCosSeries(true, ssig1, csig1, C2a);
        J12 = m0x * sig12 + (A1 * B1 - A2 * B2);
      }
    } else if ((outmask & (gReducedLength | gGeodesicScale)) != 0) {
      for (int l = 1; l <= _nC2; ++l) {
        C2a[l] = A1 * C1a[l] - A2 * C2a[l];
      }
      J12 = m0x * sig12 +
          (sinCosSeries(true, ssig2, csig2, C2a) -
           sinCosSeries(true, ssig1, csig1, C2a));
    }

    if ((outmask & gReducedLength) != 0) {
      m0 = m0x;
      m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12;
    }

    if ((outmask & gGeodesicScale) != 0) {
      final csig12 = csig1 * csig2 + ssig1 * ssig2;
      final t =
          _ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
      M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
      M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
    }

    return (s12b: s12b, m12b: m12b, m0: m0, M12: M12, M21: M21);
  }

  // ── InverseStart ──────────────────────────────────────────────────────────

  ({
    double sig12,
    double salp1, double calp1,
    double salp2, double calp2,
    double dnm
  }) _inverseStart(
      double sbet1, double cbet1, double dn1,
      double sbet2, double cbet2, double dn2,
      double lam12, double slam12, double clam12,
      List<double> C1a, List<double> C2a) {
    double sig12 = -1.0;
    double salp1 = 0.0, calp1 = 0.0;
    double salp2 = double.nan, calp2 = double.nan;
    double dnm = double.nan;

    final sbet12 = sbet2 * cbet1 - cbet2 * sbet1;
    final cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
    final sbet12a = sbet2 * cbet1 + cbet2 * sbet1;

    final shortline = cbet12 >= 0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
    double somg12, comg12;
    if (shortline) {
      double sbetm2 = GeoMath.sq(sbet1 + sbet2);
      sbetm2 /= sbetm2 + GeoMath.sq(cbet1 + cbet2);
      dnm = math.sqrt(1.0 + _ep2 * sbetm2);
      final omg12 = lam12 / (_f1 * dnm);
      somg12 = math.sin(omg12);
      comg12 = math.cos(omg12);
    } else {
      somg12 = slam12;
      comg12 = clam12;
    }

    salp1 = cbet2 * somg12;
    calp1 = comg12 >= 0
        ? sbet12 + cbet2 * sbet1 * GeoMath.sq(somg12) / (1.0 + comg12)
        : sbet12a - cbet2 * sbet1 * GeoMath.sq(somg12) / (1.0 - comg12);

    final ssig12 = GeoMath.hypot(salp1, calp1);
    final csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

    if (shortline && ssig12 < _etol2) {
      salp2 = cbet1 * somg12;
      calp2 = sbet12 - cbet1 * sbet2 *
          (comg12 >= 0
              ? GeoMath.sq(somg12) / (1.0 + comg12)
              : 1.0 - comg12);
      final t = GeoMath.hypot(salp2, calp2);
      salp2 /= t; calp2 /= t;
      sig12 = math.atan2(ssig12, csig12);
    } else if (_n.abs() > 0.1 ||
               csig12 >= 0.0 ||
               ssig12 >= 6.0 * _n.abs() * math.pi * GeoMath.sq(cbet1)) {
      // zeroth order spherical approximation OK — nothing to do
    } else {
      final lam12x = math.atan2(-slam12, -clam12);
      double x, y, lamscale, betscale;
      if (f >= 0.0) {
        final k2 = GeoMath.sq(sbet1) * _ep2;
        final eps = k2 / (2.0 * (1.0 + math.sqrt(1.0 + k2)) + k2);
        lamscale = f * cbet1 * _A3f(eps) * math.pi;
        betscale = lamscale * cbet1;
        x = lam12x / lamscale;
        y = sbet12a / betscale;
      } else {
        final cbet12a = cbet2 * cbet1 - sbet2 * sbet1;
        final bet12a = math.atan2(sbet12a, cbet12a);
        final lv = _lengths(_n, math.pi + bet12a,
            sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
            cbet1, cbet2, gReducedLength, C1a, C2a);
        final m12b = lv.m12b;
        final m0v = lv.m0;
        x = -1.0 + m12b / (cbet1 * cbet2 * m0v * math.pi);
        betscale = x < -0.01
            ? sbet12a / x
            : -f * GeoMath.sq(cbet1) * math.pi;
        lamscale = betscale / cbet1;
        y = lam12 / lamscale;
      }

      if (y > -_tol1 && x > -1.0 - _xthresh) {
        if (f >= 0.0) {
          salp1 = math.min(1.0, -x);
          calp1 = -math.sqrt(1.0 - GeoMath.sq(salp1));
        } else {
          calp1 = math.max(x > -_tol1 ? 0.0 : -1.0, x);
          salp1 = math.sqrt(1.0 - GeoMath.sq(calp1));
        }
      } else {
        final k = _astroid(x, y);
        final omg12a = lamscale *
            (f >= 0.0 ? -x * k / (1.0 + k) : -y * (1.0 + k) / k);
        somg12 = math.sin(omg12a);
        comg12 = -math.cos(omg12a);
        salp1 = cbet2 * somg12;
        calp1 = sbet12a - cbet2 * sbet1 * GeoMath.sq(somg12) / (1.0 - comg12);
      }
    }

    if (!(salp1 <= 0.0)) {
      final t = GeoMath.hypot(salp1, calp1);
      salp1 /= t; calp1 /= t;
    } else {
      salp1 = 1.0; calp1 = 0.0;
    }
    return (sig12: sig12, salp1: salp1, calp1: calp1,
            salp2: salp2, calp2: calp2, dnm: dnm);
  }

  // ── Lambda12 ──────────────────────────────────────────────────────────────

  ({
    double lam12, double salp2, double calp2, double sig12,
    double ssig1, double csig1, double ssig2, double csig2,
    double eps, double domg12, double dlam12
  }) _lambda12(
      double sbet1, double cbet1, double dn1,
      double sbet2, double cbet2, double dn2,
      double salp1, double calp1, double slam120, double clam120,
      bool diffp,
      List<double> C1a, List<double> C2a, List<double> C3a) {
    if (sbet1 == 0.0 && calp1 == 0.0) calp1 = -_tiny;

    final salp0 = salp1 * cbet1;
    final calp0 = GeoMath.hypot(calp1, salp1 * sbet1);

    double ssig1 = sbet1;
    final somg1 = salp0 * sbet1;
    double csig1 = calp1 * cbet1;
    final comg1 = calp1 * cbet1;
    var t = GeoMath.hypot(ssig1, csig1);
    ssig1 /= t; csig1 /= t;

    double salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
    double calp2 = cbet2 != cbet1 || sbet2.abs() != -sbet1
        ? math.sqrt(GeoMath.sq(calp1 * cbet1) +
                (cbet1 < -sbet1
                    ? (cbet2 - cbet1) * (cbet1 + cbet2)
                    : (sbet1 - sbet2) * (sbet1 + sbet2))) /
            cbet2
        : calp1.abs();

    double ssig2 = sbet2;
    final somg2 = salp0 * sbet2;
    double csig2 = calp2 * cbet2;
    final comg2 = calp2 * cbet2;
    t = GeoMath.hypot(ssig2, csig2);
    ssig2 /= t; csig2 /= t;

    final sig12 = math.atan2(
        math.max(0.0, csig1 * ssig2 - ssig1 * csig2),
        csig1 * csig2 + ssig1 * ssig2);

    final somg12 = math.max(0.0, comg1 * somg2 - somg1 * comg2);
    final comg12 = comg1 * comg2 + somg1 * somg2;
    final eta = math.atan2(
        somg12 * clam120 - comg12 * slam120,
        comg12 * clam120 + somg12 * slam120);

    final k2 = GeoMath.sq(calp0) * _ep2;
    final eps = k2 / (2.0 * (1.0 + math.sqrt(1.0 + k2)) + k2);
    _C3f(eps, C3a);
    final B312 = sinCosSeries(true, ssig2, csig2, C3a) -
        sinCosSeries(true, ssig1, csig1, C3a);
    final domg12 = -f * _A3f(eps) * salp0 * (sig12 + B312);
    final lam12 = eta + domg12;

    double dlam12 = double.nan;
    if (diffp) {
      if (calp2 == 0.0) {
        dlam12 = -2.0 * _f1 * dn1 / sbet1;
      } else {
        final lv = _lengths(eps, sig12,
            ssig1, csig1, dn1, ssig2, csig2, dn2,
            cbet1, cbet2, gReducedLength, C1a, C2a);
        dlam12 = lv.m12b;
        dlam12 *= _f1 / (calp2 * cbet2);
      }
    }

    return (lam12: lam12, salp2: salp2, calp2: calp2, sig12: sig12,
            ssig1: ssig1, csig1: csig1, ssig2: ssig2, csig2: csig2,
            eps: eps, domg12: domg12, dlam12: dlam12);
  }

  // ── InverseInt ────────────────────────────────────────────────────────────

  _InverseResult _inverseInt(
      double lat1, double lon1, double lat2, double lon2, int outmask) {
    final vals = GeodesicData();
    vals.lat1 = lat1 = GeoMath.latFix(lat1);
    vals.lat2 = lat2 = GeoMath.latFix(lat2);
    lat1 = GeoMath.angRound(lat1);
    lat2 = GeoMath.angRound(lat2);

    final ad = GeoMath.angDiff(lon1, lon2);
    double lon12 = ad.$1;
    double lon12s = ad.$2;

    if ((outmask & gLongUnroll) != 0) {
      vals.lon1 = lon1;
      vals.lon2 = (lon1 + lon12) + lon12s;
    } else {
      vals.lon1 = GeoMath.angNormalize(lon1);
      vals.lon2 = GeoMath.angNormalize(lon2);
    }

    double lonsign = GeoMath.copysign(1.0, lon12);
    lon12 *= lonsign;
    lon12s *= lonsign;
    final lam12 = lon12 * GeoMath.degree;
    final sc = GeoMath.sincosde(lon12, lon12s);
    double slam12 = sc.$1, clam12 = sc.$2;
    lon12s = (180.0 - lon12) - lon12s;

    double swappD = (lat1.abs() < lat2.abs() || lat2.isNaN) ? -1.0 : 1.0;
    if (swappD < 0) {
      lonsign *= -1.0;
      final tmp1 = lat2; lat2 = lat1; lat1 = tmp1;
    }
    double latsign = GeoMath.copysign(1.0, -lat1);
    lat1 *= latsign;
    lat2 *= latsign;

    var sc1 = GeoMath.sincosd(lat1);
    double sbet1 = _f1 * sc1.$1, cbet1 = sc1.$2;
    var t2 = GeoMath.hypot(sbet1, cbet1);
    sbet1 /= t2; cbet1 /= t2;
    cbet1 = math.max(_tiny, cbet1);

    var sc2 = GeoMath.sincosd(lat2);
    double sbet2 = _f1 * sc2.$1, cbet2 = sc2.$2;
    t2 = GeoMath.hypot(sbet2, cbet2);
    sbet2 /= t2; cbet2 /= t2;
    cbet2 = math.max(_tiny, cbet2);

    if (cbet1 < -sbet1) {
      if (cbet2 == cbet1) sbet2 = GeoMath.copysign(sbet1, sbet2);
    } else {
      if (sbet2.abs() == -sbet1) cbet2 = cbet1;
    }

    final dn1 = math.sqrt(1.0 + _ep2 * GeoMath.sq(sbet1));
    final dn2 = math.sqrt(1.0 + _ep2 * GeoMath.sq(sbet2));

    final C1a = List<double>.filled(_nC1 + 1, 0.0);
    final C2a = List<double>.filled(_nC2 + 1, 0.0);
    final C3a = List<double>.filled(_nC3, 0.0);

    bool meridian = lat1 == -90.0 || slam12 == 0.0;
    double s12x = 0.0, m12x = 0.0;
    double sig12 = 0.0, ssig1 = 0.0, csig1 = 0.0;
    double ssig2 = 0.0, csig2 = 0.0;
    double salp1 = 0.0, calp1 = 0.0, salp2 = 0.0, calp2 = 0.0;
    double eps = 0.0, omg12 = 0.0, domg12 = 0.0;

    if (meridian) {
      calp1 = clam12; salp1 = slam12;
      calp2 = 1.0; salp2 = 0.0;
      ssig1 = sbet1; csig1 = calp1 * cbet1;
      ssig2 = sbet2; csig2 = calp2 * cbet2;

      sig12 = math.atan2(
          math.max(0.0, csig1 * ssig2 - ssig1 * csig2),
          csig1 * csig2 + ssig1 * ssig2);
      final lv = _lengths(_n, sig12,
          ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
          outmask | gDistance | gReducedLength, C1a, C2a);
      s12x = lv.s12b;
      m12x = lv.m12b;
      if ((outmask & gGeodesicScale) != 0) {
        vals.M12 = lv.M12;
        vals.M21 = lv.M21;
      }
      if (sig12 < _tol2 || m12x >= 0.0) {
        if (sig12 < 3.0 * _tiny ||
            (sig12 < _tol0 && (s12x < 0.0 || m12x < 0.0))) {
          sig12 = m12x = s12x = 0.0;
        }
        m12x *= _b;
        s12x *= _b;
        vals.a12 = sig12 / GeoMath.degree;
      } else {
        meridian = false;
      }
    }

    double somg12 = 2.0;
    if (!meridian &&
        sbet1 == 0.0 &&
        (f <= 0.0 || lon12s >= f * 180.0)) {
      calp1 = calp2 = 0.0;
      salp1 = salp2 = 1.0;
      s12x = a * lam12;
      sig12 = omg12 = lam12 / _f1;
      m12x = _b * math.sin(sig12);
      if ((outmask & gGeodesicScale) != 0) {
        vals.M12 = vals.M21 = math.cos(sig12);
      }
      vals.a12 = lon12 / _f1;
    } else if (!meridian) {
      final iv = _inverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
          lam12, slam12, clam12, C1a, C2a);
      sig12 = iv.sig12;
      salp1 = iv.salp1; calp1 = iv.calp1;

      if (sig12 >= 0.0) {
        salp2 = iv.salp2; calp2 = iv.calp2;
        final dnm = iv.dnm;
        s12x = sig12 * _b * dnm;
        m12x = GeoMath.sq(dnm) * _b * math.sin(sig12 / dnm);
        if ((outmask & gGeodesicScale) != 0) {
          vals.M12 = vals.M21 = math.cos(sig12 / dnm);
        }
        vals.a12 = sig12 / GeoMath.degree;
        omg12 = lam12 / (_f1 * dnm);
      } else {
        int numit = 0;
        double salp1a = _tiny, calp1a = 1.0;
        double salp1b = _tiny, calp1b = -1.0;
        bool tripn = false, tripb = false;

        while (true) {
          final lv = _lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
              salp1, calp1, slam12, clam12, numit < _maxit1,
              C1a, C2a, C3a);
          double v = lv.lam12;
          salp2 = lv.salp2; calp2 = lv.calp2;
          sig12 = lv.sig12;
          ssig1 = lv.ssig1; csig1 = lv.csig1;
          ssig2 = lv.ssig2; csig2 = lv.csig2;
          eps = lv.eps;
          domg12 = lv.domg12;
          final dv = lv.dlam12;

          if (tripb ||
              !(v.abs() >= (tripn ? 8.0 : 1.0) * _tol0) ||
              numit == _maxit2) break;

          if (v > 0.0 && (numit < _maxit1 || calp1 / salp1 > calp1b / salp1b)) {
            salp1b = salp1; calp1b = calp1;
          } else if (v < 0.0 &&
              (numit < _maxit1 || calp1 / salp1 < calp1a / salp1a)) {
            salp1a = salp1; calp1a = calp1;
          }

          if (numit < _maxit1 && dv > 0.0) {
            final dalp1 = -v / dv;
            if (dalp1.abs() < math.pi) {
              final sdalp1 = math.sin(dalp1);
              final cdalp1 = math.cos(dalp1);
              final nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
              if (nsalp1 > 0.0) {
                calp1 = calp1 * cdalp1 - salp1 * sdalp1;
                salp1 = nsalp1;
                final tn = GeoMath.hypot(salp1, calp1);
                salp1 /= tn; calp1 /= tn;
                tripn = v.abs() <= 16.0 * _tol0;
                ++numit;
                continue;
              }
            }
          }
          salp1 = (salp1a + salp1b) / 2.0;
          calp1 = (calp1a + calp1b) / 2.0;
          final tn = GeoMath.hypot(salp1, calp1);
          salp1 /= tn; calp1 /= tn;
          tripn = false;
          tripb = (salp1a - salp1).abs() + (calp1a - calp1) < _tolb ||
                  (salp1 - salp1b).abs() + (calp1 - calp1b) < _tolb;
          ++numit;
        }

        final lengthmask = outmask |
            ((outmask & (gReducedLength | gGeodesicScale)) != 0
                ? gDistance
                : gNone);
        final lv = _lengths(eps, sig12,
            ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
            lengthmask, C1a, C2a);
        s12x = lv.s12b;
        m12x = lv.m12b;
        if ((outmask & gGeodesicScale) != 0) {
          vals.M12 = lv.M12;
          vals.M21 = lv.M21;
        }
        m12x *= _b;
        s12x *= _b;
        vals.a12 = sig12 / GeoMath.degree;
        if ((outmask & gArea) != 0) {
          final sdomg12 = math.sin(domg12);
          final cdomg12 = math.cos(domg12);
          somg12 = slam12 * cdomg12 - clam12 * sdomg12;
          final comg12t = clam12 * cdomg12 + slam12 * sdomg12;
          omg12 = math.atan2(somg12, comg12t);
        }
      }
    }

    if ((outmask & gDistance) != 0) {
      vals.s12 = 0.0 + s12x;
    }
    if ((outmask & gReducedLength) != 0) {
      vals.m12 = 0.0 + m12x;
    }

    if ((outmask & gArea) != 0) {
      final salp0 = salp1 * cbet1;
      final calp0 = GeoMath.hypot(calp1, salp1 * sbet1);
      double areaS12 = 0.0;
      if (calp0 != 0.0 && salp0 != 0.0) {
        double ssig1t = sbet1, csig1t = calp1 * cbet1;
        double ssig2t = sbet2, csig2t = calp2 * cbet2;
        final k2 = GeoMath.sq(calp0) * _ep2;
        final epst = k2 / (2.0 * (1.0 + math.sqrt(1.0 + k2)) + k2);
        final A4 = GeoMath.sq(a) * calp0 * salp0 * _e2;
        var tn = GeoMath.hypot(ssig1t, csig1t);
        ssig1t /= tn; csig1t /= tn;
        tn = GeoMath.hypot(ssig2t, csig2t);
        ssig2t /= tn; csig2t /= tn;
        final C4a = List<double>.filled(_nC4, 0.0);
        _C4f(epst, C4a);
        final B41 = sinCosSeries(false, ssig1t, csig1t, C4a);
        final B42 = sinCosSeries(false, ssig2t, csig2t, C4a);
        areaS12 = A4 * (B42 - B41);
      }
      if (!meridian && somg12 == 2.0) {
        somg12 = math.sin(omg12);
      }

      final comg12 = math.cos(omg12);
      double alp12;
      if (!meridian && comg12 > -0.7071 && sbet2 - sbet1 < 1.75) {
        final domg12t = 1.0 + comg12;
        final dbet1 = 1.0 + cbet1;
        final dbet2 = 1.0 + cbet2;
        alp12 = 2.0 * math.atan2(
            somg12 * (sbet1 * dbet2 + sbet2 * dbet1),
            domg12t * (sbet1 * sbet2 + dbet1 * dbet2));
      } else {
        double salp12 = salp2 * calp1 - calp2 * salp1;
        double calp12 = calp2 * calp1 + salp2 * salp1;
        if (salp12 == 0.0 && calp12 < 0.0) {
          salp12 = _tiny * calp1;
          calp12 = -1.0;
        }
        alp12 = math.atan2(salp12, calp12);
      }
      areaS12 += _c2 * alp12;
      areaS12 *= swappD * lonsign * latsign;
      areaS12 += 0.0;
      vals.S12 = areaS12;
    }

    if (swappD < 0.0) {
      final ts1 = salp2; salp2 = salp1; salp1 = ts1;
      final tc1 = calp2; calp2 = calp1; calp1 = tc1;
      if ((outmask & gGeodesicScale) != 0) {
        final tm = vals.M21; vals.M21 = vals.M12; vals.M12 = tm;
      }
    }

    salp1 *= swappD * lonsign; calp1 *= swappD * latsign;
    salp2 *= swappD * lonsign; calp2 *= swappD * latsign;

    return _InverseResult(vals, salp1, calp1, salp2, calp2);
  }

  // ── Public API ────────────────────────────────────────────────────────────

  /// Solve the geodesic inverse problem.
  GeodesicData inverse(double lat1, double lon1, double lat2, double lon2,
      [int outmask = gStandard]) {
    if (outmask == gLongUnroll) outmask |= gStandard;
    outmask &= gOutMask;
    final r = _inverseInt(lat1, lon1, lat2, lon2, outmask);
    final vals = r.vals;
    if ((outmask & gAzimuth) != 0) {
      vals.azi1 = GeoMath.atan2d(r.salp1, r.calp1);
      vals.azi2 = GeoMath.atan2d(r.salp2, r.calp2);
    }
    return vals;
  }

  /// Solve the general direct geodesic problem.
  GeodesicData genDirect(double lat1, double lon1, double azi1,
      bool arcmode, double s12_a12,
      [int outmask = gStandard]) {
    if (outmask == gLongUnroll) outmask |= gStandard;
    if (!arcmode) outmask |= gDistanceIn;
    final line = GeodesicLine(this, lat1, lon1, azi1, outmask);
    return line.genPosition(arcmode, s12_a12, outmask);
  }

  /// Solve the direct geodesic problem (distance in meters).
  GeodesicData direct(double lat1, double lon1, double azi1, double s12,
      [int outmask = gStandard]) {
    return genDirect(lat1, lon1, azi1, false, s12, outmask);
  }

  /// Solve the direct geodesic problem with arc length in degrees.
  GeodesicData arcDirect(double lat1, double lon1, double azi1, double a12,
      [int outmask = gStandard]) {
    return genDirect(lat1, lon1, azi1, true, a12, outmask);
  }

  /// Create a GeodesicLine.
  GeodesicLine line(double lat1, double lon1, double azi1,
      [int caps = gStandard | gDistanceIn]) {
    return GeodesicLine(this, lat1, lon1, azi1, caps);
  }

  /// Define a GeodesicLine via the direct problem (distance).
  GeodesicLine directLine(double lat1, double lon1, double azi1, double s12,
      [int caps = gStandard | gDistanceIn]) {
    return _genDirectLine(lat1, lon1, azi1, false, s12, caps);
  }

  /// Define a GeodesicLine via the direct problem (arc length).
  GeodesicLine arcDirectLine(double lat1, double lon1, double azi1, double a12,
      [int caps = gStandard | gDistanceIn]) {
    return _genDirectLine(lat1, lon1, azi1, true, a12, caps);
  }

  GeodesicLine _genDirectLine(double lat1, double lon1, double azi1,
      bool arcmode, double s12_a12,
      [int caps = gStandard | gDistanceIn]) {
    if (!arcmode) caps |= gDistanceIn;
    final t = GeodesicLine(this, lat1, lon1, azi1, caps);
    t.genSetDistance(arcmode, s12_a12);
    return t;
  }

  /// Define a GeodesicLine via the inverse problem.
  GeodesicLine inverseLine(double lat1, double lon1, double lat2, double lon2,
      [int caps = gStandard | gDistanceIn]) {
    final r = _inverseInt(lat1, lon1, lat2, lon2, gArc);
    final azi1 = GeoMath.atan2d(r.salp1, r.calp1);
    int capst = caps;
    if ((capst & (gOutMask & 0x7F80 & gDistanceIn)) != 0) capst |= gDistance;
    final t = GeodesicLine(this, lat1, lon1, azi1, capst, r.salp1, r.calp1);
    t.setArc(r.vals.a12);
    return t;
  }

  /// Create a PolygonArea for this ellipsoid.
  PolygonArea polygon([bool polyline = false]) => PolygonArea(this, polyline);
}

// ══════════════════════════════════════════════════════════════════════════════
// GeodesicLine (defined in same file to avoid circular import)
// ══════════════════════════════════════════════════════════════════════════════

/// Performs geodesic calculations along a given geodesic line.
class GeodesicLine {
  final double a;
  final double f;
  final double _b;
  final double _c2;
  final double _f1;
  int caps;

  late double lat1, lon1, azi1;
  late double salp1, calp1;
  double s13 = double.nan;
  double a13 = double.nan;

  late double _dn1, _salp0, _calp0;
  late double _ssig1, _csig1, _somg1, _comg1, _k2;

  late double _A1m1, _B11, _stau1, _ctau1;
  late List<double> _C1a, _C1pa, _C2a, _C3a, _C4a;
  late double _A2m1, _B21, _A3c, _B31, _A4, _B41;

  GeodesicLine(Geodesic geod, double lat1, double lon1, double azi1,
      [int caps = gStandard | gDistanceIn,
      double? salp1in, double? calp1in])
      : a = geod.a,
        f = geod.f,
        _b = geod._b,
        _c2 = geod._c2,
        _f1 = geod._f1,
        caps = (caps == 0 ? gStandard | gDistanceIn : caps) |
            gLatitude | gAzimuth | gLongUnroll {
    this.lat1 = GeoMath.latFix(lat1);
    this.lon1 = lon1;
    if (salp1in == null || calp1in == null) {
      this.azi1 = GeoMath.angNormalize(azi1);
      final t = GeoMath.sincosd(GeoMath.angRound(this.azi1));
      this.salp1 = t.$1; this.calp1 = t.$2;
    } else {
      this.azi1 = azi1;
      this.salp1 = salp1in; this.calp1 = calp1in;
    }

    final t1 = GeoMath.sincosd(GeoMath.angRound(this.lat1));
    double sbet1 = _f1 * t1.$1, cbet1 = t1.$2;
    var tn = GeoMath.hypot(sbet1, cbet1);
    sbet1 /= tn; cbet1 /= tn;
    cbet1 = math.max(_tiny, cbet1);
    _dn1 = math.sqrt(1.0 + geod._ep2 * GeoMath.sq(sbet1));

    _salp0 = this.salp1 * cbet1;
    _calp0 = GeoMath.hypot(this.calp1, this.salp1 * sbet1);

    _ssig1 = sbet1;
    _somg1 = _salp0 * sbet1;
    _csig1 = _comg1 =
        (sbet1 != 0.0 || this.calp1 != 0.0) ? cbet1 * this.calp1 : 1.0;
    tn = GeoMath.hypot(_ssig1, _csig1);
    _ssig1 /= tn; _csig1 /= tn;

    _k2 = GeoMath.sq(_calp0) * geod._ep2;
    final eps = _k2 / (2.0 * (1.0 + math.sqrt(1.0 + _k2)) + _k2);

    if ((this.caps & capC1) != 0) {
      _A1m1 = _A1m1f(eps);
      _C1a = List<double>.filled(_nC1 + 1, 0.0);
      _C1f(eps, _C1a);
      _B11 = sinCosSeries(true, _ssig1, _csig1, _C1a);
      final s = math.sin(_B11), c = math.cos(_B11);
      _stau1 = _ssig1 * c + _csig1 * s;
      _ctau1 = _csig1 * c - _ssig1 * s;
    }

    if ((this.caps & capC1p) != 0) {
      _C1pa = List<double>.filled(_nC1p + 1, 0.0);
      _C1pf(eps, _C1pa);
    }

    if ((this.caps & capC2) != 0) {
      _A2m1 = _A2m1f(eps);
      _C2a = List<double>.filled(_nC2 + 1, 0.0);
      _C2f(eps, _C2a);
      _B21 = sinCosSeries(true, _ssig1, _csig1, _C2a);
    }

    if ((this.caps & capC3) != 0) {
      _C3a = List<double>.filled(_nC3, 0.0);
      geod._C3f(eps, _C3a);
      _A3c = -f * _salp0 * geod._A3f(eps);
      _B31 = sinCosSeries(true, _ssig1, _csig1, _C3a);
    }

    if ((this.caps & capC4) != 0) {
      _C4a = List<double>.filled(_nC4, 0.0);
      geod._C4f(eps, _C4a);
      _A4 = GeoMath.sq(a) * _calp0 * _salp0 * geod._e2;
      _B41 = sinCosSeries(false, _ssig1, _csig1, _C4a);
    }
  }

  /// Compute position (general: arc or distance mode).
  GeodesicData genPosition(bool arcmode, double s12_a12, [int outmask = gStandard]) {
    if (outmask == gLongUnroll) outmask |= gStandard;
    outmask &= caps & gOutMask;
    final vals = GeodesicData();
    vals.lat1 = lat1;
    vals.azi1 = azi1;
    vals.lon1 = (outmask & gLongUnroll) != 0
        ? lon1
        : GeoMath.angNormalize(lon1);
    if (arcmode) {
      vals.a12 = s12_a12;
    } else {
      vals.s12 = s12_a12;
    }
    if (!(arcmode || ((caps & gDistanceIn & gOutMask) != 0))) {
      vals.a12 = double.nan;
      return vals;
    }

    double B12 = 0.0, AB1 = 0.0;
    double sig12 = 0.0, ssig12 = 0.0, csig12 = 1.0;

    if (arcmode) {
      sig12 = s12_a12 * GeoMath.degree;
      final sc = GeoMath.sincosd(s12_a12);
      ssig12 = sc.$1; csig12 = sc.$2;
    } else {
      final tau12 = s12_a12 / (_b * (1.0 + _A1m1));
      final s = math.sin(tau12), c = math.cos(tau12);
      B12 = -sinCosSeries(true,
          _stau1 * c + _ctau1 * s,
          _ctau1 * c - _stau1 * s,
          _C1pa);
      sig12 = tau12 - (B12 - _B11);
      ssig12 = math.sin(sig12); csig12 = math.cos(sig12);
      if (f.abs() > 0.01) {
        double ssig2t = _ssig1 * csig12 + _csig1 * ssig12;
        double csig2t = _csig1 * csig12 - _ssig1 * ssig12;
        B12 = sinCosSeries(true, ssig2t, csig2t, _C1a);
        final serr = (1.0 + _A1m1) * (sig12 + (B12 - _B11)) - s12_a12 / _b;
        sig12 = sig12 - serr / math.sqrt(1.0 + _k2 * GeoMath.sq(ssig2t));
        ssig12 = math.sin(sig12); csig12 = math.cos(sig12);
      }
    }

    double ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
    double csig2 = _csig1 * csig12 - _ssig1 * ssig12;
    final dn2 = math.sqrt(1.0 + _k2 * GeoMath.sq(ssig2));

    if ((outmask & (gDistance | gReducedLength | gGeodesicScale)) != 0) {
      if (arcmode || f.abs() > 0.01) {
        B12 = sinCosSeries(true, ssig2, csig2, _C1a);
      }
      AB1 = (1.0 + _A1m1) * (B12 - _B11);
    }

    final sbet2 = _calp0 * ssig2;
    double cbet2 = GeoMath.hypot(_salp0, _calp0 * csig2);
    if (cbet2 == 0.0) cbet2 = csig2 = _tiny;
    final salp2 = _salp0, calp2 = _calp0 * csig2;

    if (arcmode && (outmask & gDistance) != 0) {
      vals.s12 = _b * ((1.0 + _A1m1) * sig12 + AB1);
    }

    if ((outmask & gLongitude) != 0) {
      final somg2 = _salp0 * ssig2, comg2 = csig2;
      final E = GeoMath.copysign(1.0, _salp0);
      double omg12;
      if ((outmask & gLongUnroll) != 0) {
        omg12 = E * (sig12 -
            (math.atan2(ssig2, csig2) - math.atan2(_ssig1, _csig1)) +
            (math.atan2(E * somg2, comg2) - math.atan2(E * _somg1, _comg1)));
      } else {
        omg12 = math.atan2(
            somg2 * _comg1 - comg2 * _somg1,
            comg2 * _comg1 + somg2 * _somg1);
      }
      final lam12 = omg12 + _A3c * (sig12 +
          (sinCosSeries(true, ssig2, csig2, _C3a) - _B31));
      final lon12 = lam12 / GeoMath.degree;
      vals.lon2 = (outmask & gLongUnroll) != 0
          ? lon1 + lon12
          : GeoMath.angNormalize(GeoMath.angNormalize(lon1) +
              GeoMath.angNormalize(lon12));
    }

    if ((outmask & gLatitude) != 0) {
      vals.lat2 = GeoMath.atan2d(sbet2, _f1 * cbet2);
    }

    if ((outmask & gAzimuth) != 0) {
      vals.azi2 = GeoMath.atan2d(salp2, calp2);
    }

    if ((outmask & (gReducedLength | gGeodesicScale)) != 0) {
      final B22 = sinCosSeries(true, ssig2, csig2, _C2a);
      final AB2 = (1.0 + _A2m1) * (B22 - _B21);
      final J12 = (_A1m1 - _A2m1) * sig12 + (AB1 - AB2);
      if ((outmask & gReducedLength) != 0) {
        vals.m12 = _b *
            ((dn2 * (_csig1 * ssig2) - _dn1 * (_ssig1 * csig2)) -
             _csig1 * csig2 * J12);
      }
      if ((outmask & gGeodesicScale) != 0) {
        final t = _k2 * (ssig2 - _ssig1) * (ssig2 + _ssig1) / (_dn1 + dn2);
        vals.M12 = csig12 + (t * ssig2 - csig2 * J12) * _ssig1 / _dn1;
        vals.M21 = csig12 - (t * _ssig1 - _csig1 * J12) * ssig2 / dn2;
      }
    }

    if ((outmask & gArea) != 0) {
      final B42 = sinCosSeries(false, ssig2, csig2, _C4a);
      double salp12, calp12;
      if (_calp0 == 0.0 || _salp0 == 0.0) {
        salp12 = salp2 * calp1 - calp2 * salp1;
        calp12 = calp2 * calp1 + salp2 * salp1;
      } else {
        salp12 = _calp0 * _salp0 *
            (csig12 <= 0.0
                ? _csig1 * (1.0 - csig12) + ssig12 * _ssig1
                : ssig12 * (_csig1 * ssig12 / (1.0 + csig12) + _ssig1));
        calp12 = GeoMath.sq(_salp0) + GeoMath.sq(_calp0) * _csig1 * csig2;
      }
      vals.S12 = _c2 * math.atan2(salp12, calp12) + _A4 * (B42 - _B41);
    }

    if (!arcmode) vals.a12 = sig12 / GeoMath.degree;
    return vals;
  }

  /// Position at distance [s12] meters.
  GeodesicData position(double s12, [int outmask = gStandard]) {
    return genPosition(false, s12, outmask);
  }

  /// Position at arc length [a12] degrees.
  GeodesicData arcPosition(double a12, [int outmask = gStandard]) {
    return genPosition(true, a12, outmask);
  }

  void genSetDistance(bool arcmode, double s13_a13) {
    if (arcmode) setArc(s13_a13); else setDistance(s13_a13);
  }

  void setDistance(double s13) {
    this.s13 = s13;
    final r = genPosition(false, s13, gArc);
    a13 = 0.0 + r.a12;
  }

  void setArc(double a13) {
    this.a13 = a13;
    final r = genPosition(true, a13, gDistance);
    s13 = 0.0 + r.s12;
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// PolygonArea
// ══════════════════════════════════════════════════════════════════════════════

int _transit(double lon1, double lon2) {
  final lon12 = GeoMath.angDiff(lon1, lon2).$1;
  final l1 = GeoMath.angNormalize(lon1);
  final l2 = GeoMath.angNormalize(lon2);
  return lon12 > 0.0 &&
          ((l1 < 0.0 && l2 >= 0.0) || (l1 > 0.0 && l2 == 0.0))
      ? 1
      : (lon12 < 0.0 && l1 >= 0.0 && l2 < 0.0 ? -1 : 0);
}

int _transitDirect(double lon1, double lon2) {
  lon1 = lon1 % 720.0;
  lon2 = lon2 % 720.0;
  return ((lon2 >= 0.0 && lon2 < 360.0) || lon2 < -360.0 ? 0 : 1) -
      ((lon1 >= 0.0 && lon1 < 360.0) || lon1 < -360.0 ? 0 : 1);
}

double _areaReduceA(
    Accumulator area, double area0, int crossings, bool reverse, bool sign) {
  area.remainder(area0);
  if ((crossings & 1) != 0) {
    area.add((area.sum() < 0.0 ? 1.0 : -1.0) * area0 / 2.0);
  }
  if (!reverse) area.negate();
  if (sign) {
    if (area.sum() > area0 / 2.0) area.add(-area0);
    else if (area.sum() <= -area0 / 2.0) area.add(area0);
  } else {
    if (area.sum() >= area0) area.add(-area0);
    else if (area.sum() < 0.0) area.add(area0);
  }
  return 0.0 + area.sum();
}

double _areaReduceB(
    double area, double area0, int crossings, bool reverse, bool sign) {
  area = GeoMath.remainder(area, area0);
  if ((crossings & 1) != 0) {
    area += (area < 0.0 ? 1.0 : -1.0) * area0 / 2.0;
  }
  if (!reverse) area *= -1.0;
  if (sign) {
    if (area > area0 / 2.0) area -= area0;
    else if (area <= -area0 / 2.0) area += area0;
  } else {
    if (area >= area0) area -= area0;
    else if (area < 0.0) area += area0;
  }
  return 0.0 + area;
}

/// Result of PolygonArea.compute.
class PolygonResult {
  final int number;
  final double perimeter;
  final double area;
  PolygonResult(this.number, this.perimeter, this.area);
}

/// Computes geodesic polygon area and perimeter.
class PolygonArea {
  final Geodesic _geod;
  final double a;
  final double f;
  final double _area0;
  final bool polyline;
  late int _mask;

  Accumulator? _areasum;
  late Accumulator _perimetersum;

  int num = 0;
  int _crossings = 0;
  double lat = double.nan, lon = double.nan;
  double _lat0 = double.nan, _lon0 = double.nan;

  PolygonArea(this._geod, [this.polyline = false])
      : a = _geod.a,
        f = _geod.f,
        _area0 = 4.0 * math.pi * _geod._c2 {
    _mask = gLatitude | gLongitude | gDistance |
        (polyline ? gNone : gArea | gLongUnroll);
    if (!polyline) _areasum = Accumulator(0.0);
    _perimetersum = Accumulator(0.0);
    clear();
  }

  void clear() {
    num = 0;
    _crossings = 0;
    if (!polyline) _areasum!.set(0.0);
    _perimetersum.set(0.0);
    _lat0 = _lon0 = lat = lon = double.nan;
  }

  void addPoint(double lat, double lon) {
    if (num == 0) {
      _lat0 = this.lat = lat;
      _lon0 = this.lon = lon;
    } else {
      final t = _geod.inverse(this.lat, this.lon, lat, lon, _mask);
      _perimetersum.add(t.s12);
      if (!polyline) {
        _areasum!.add(t.S12);
        _crossings += _transit(this.lon, lon);
      }
      this.lat = lat;
      this.lon = lon;
    }
    ++num;
  }

  void addEdge(double azi, double s) {
    if (num > 0) {
      final t = _geod.direct(lat, lon, azi, s, _mask);
      _perimetersum.add(s);
      if (!polyline) {
        _areasum!.add(t.S12);
        _crossings += _transitDirect(lon, t.lon2);
      }
      lat = t.lat2;
      lon = t.lon2;
    }
    ++num;
  }

  PolygonResult compute({bool reverse = false, bool sign = false}) {
    if (num < 2) {
      return PolygonResult(num, 0.0, polyline ? double.nan : 0.0);
    }
    if (polyline) {
      return PolygonResult(num, _perimetersum.sum(), double.nan);
    }
    final t = _geod.inverse(lat, lon, _lat0, _lon0, _mask);
    final perimeter = _perimetersum.sum(t.s12);
    final tempsum = Accumulator(_areasum!);
    tempsum.add(t.S12);
    final area = _areaReduceA(
        tempsum, _area0, _crossings + _transit(lon, _lon0), reverse, sign);
    return PolygonResult(num, perimeter, area);
  }

  PolygonResult testPoint(double lat, double lon,
      {bool reverse = false, bool sign = false}) {
    if (num == 0) {
      return PolygonResult(1, 0.0, polyline ? double.nan : 0.0);
    }
    double perimeter = _perimetersum.sum();
    double tempsum = polyline ? 0.0 : _areasum!.sum();
    int crossings = _crossings;
    for (int i = 0; i < (polyline ? 1 : 2); ++i) {
      final t = _geod.inverse(
          i == 0 ? this.lat : lat, i == 0 ? this.lon : lon,
          i != 0 ? _lat0 : lat, i != 0 ? _lon0 : lon,
          _mask);
      perimeter += t.s12;
      if (!polyline) {
        tempsum += t.S12;
        crossings += _transit(
            i == 0 ? this.lon : lon,
            i != 0 ? _lon0 : lon);
      }
    }
    if (polyline) return PolygonResult(num + 1, perimeter, double.nan);
    final area = _areaReduceB(tempsum, _area0, crossings, reverse, sign);
    return PolygonResult(num + 1, perimeter, area);
  }

  PolygonResult testEdge(double azi, double s,
      {bool reverse = false, bool sign = false}) {
    if (num == 0) return PolygonResult(0, double.nan, double.nan);
    double perimeter = _perimetersum.sum() + s;
    if (polyline) return PolygonResult(num + 1, perimeter, double.nan);

    double tempsum = _areasum!.sum();
    int crossings = _crossings;
    final t = _geod.direct(lat, lon, azi, s, _mask);
    tempsum += t.S12;
    crossings += _transitDirect(lon, t.lon2);
    crossings += _transit(t.lon2, _lon0);
    final t2 = _geod.inverse(t.lat2, t.lon2, _lat0, _lon0, _mask);
    perimeter += t2.s12;
    tempsum += t2.S12;

    final area = _areaReduceB(tempsum, _area0, crossings, reverse, sign);
    return PolygonResult(num + 1, perimeter, area);
  }
}
