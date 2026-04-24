// geodesic_line3.dart
// Dart port of GeodesicLine3.js / GeodesicLine3.{hpp,cpp}
// Part of GeographicLib (C) Charles Karney, triaxial geodesic.
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
import 'angle.dart';
import 'elliptic_function3.dart';
import 'trigfun.dart';
import 'ellipsoid3.dart';

// ignore_for_file: non_constant_identifier_names

// ---------------------------------------------------------------------------
// Module constants and helpers
// ---------------------------------------------------------------------------

const double _pi = math.pi;
const double _machEps = 2.220446049250313e-16;
const int _maxit = 300;

// BigValue: a large constant used as "infinity" for lambertian bounds
double _bigValue() => -3.0 * math.log(_machEps); // ≈ 104.2

double _sq(double x) => x * x;

bool _signbit(double x) =>
    x.isNegative || (x == 0.0 && (1.0 / x).isNegative);

double _copysign(double x, double y) =>
    x.abs() * (_signbit(y) ? -1.0 : 1.0);

double _clamp(double x, double lo, double hi) =>
    math.max(lo, math.min(hi, x));

double _bigclamp(double x, [double mult = 1.0]) {
  final z = mult * _bigValue();
  return _clamp(x, -z, z);
}

// lam(x, mult=1) = |x|>=pi/2 ? copysign(BigValue,x) : asinh(mult*tan(x))
double _lam(double x, [double mult = 1.0]) {
  return x.abs() >= _pi / 2.0
      ? _copysign(_bigValue(), x)
      : math.log(mult * math.tan(x) +
            math.sqrt(1.0 + _sq(mult * math.tan(x))));
}

// gd(x, mult=1) = atan(sinh3(x)/mult)
double _gd(double x, [double mult = 1.0]) {
  return math.atan(sinh3(x) / mult);
}

// anglam(u, mult=1): Angle(sinh3(u), mult) normalised
Angle _anglam(double u, [double mult = 1.0]) {
  final sh = sinh3(u);
  final h = math.sqrt(sh * sh + mult * mult);
  return Angle(sh / h, mult / h, 0, true);
}

// lamang(ang, mult=1) = bigclamp(asinh(mult * ang.t()))
double _lamang(Angle ang, [double mult = 1.0]) {
  if (ang.cx == 0.0) {
    return _bigclamp(_copysign(_bigValue(), mult * ang.sx));
  }
  final t = mult * ang.sx / ang.cx;
  // asinh(t): for t<0 use -asinh(-t) to avoid catastrophic cancellation
  final at = t.abs();
  final r = math.log(at + math.sqrt(1.0 + at * at));
  return _bigclamp(_copysign(r, t));
}

// mcosh(u, mult=1) = cosh3(u) if mult==1, else hypot(sinh3(u),mult)/mult
double _mcosh(double u, double mult) {
  return mult == 1.0
      ? cosh3(u)
      : math.sqrt(_sq(sinh3(u)) + _sq(mult)) / mult;
}

// modang(x, m) = atan(m * tan(x)) via Angle
double _modang(double x, double m) {
  return Angle.fromRadians(x).modang(m).radians();
}

// remxAng: reduce ang to [-pi/2, pi/2), returns (val, n)
({double val, int n}) _remxAng(Angle ang, [bool alt = false]) {
  final x = ang.copy();
  var m = 0;
  if (_signbit(x.cx)) {
    x.subEq(Angle.cardinal(2));
    ++m;
  }
  if (x.cx == 0.0) {
    if (alt && _signbit(x.sx)) {
      x.addEq(Angle.cardinal(2));
      --m;
    } else if (!alt && !_signbit(x.sx)) {
      x.subEq(Angle.cardinal(2));
      ++m;
    }
  }
  return (val: x.radians0(), n: m + 2 * x.nx);
}

// remxReal: remainder of x mod y in [-y/2, y/2)
({double val, int n}) _remxReal(double x, double y, [bool alt = false]) {
  var n = (x / y).roundToDouble().toInt();
  var z = x - n * y;
  if (alt) {
    if (z == -y / 2.0) z = y / 2.0;
  } else {
    if (z == y / 2.0) z = -y / 2.0;
  }
  return (val: z, n: n);
}

// biaxspecial: special biaxial limit case?
bool _biaxspecial(Geodesic3Ref tg, double gamma) {
  return (tg.t.k2 == 0.0 || tg.t.kp2 == 0.0) &&
      gamma != 0.0 &&
      gamma.abs() < tg.ellipthresh;
}

// ---------------------------------------------------------------------------
// Forward declaration for Geodesic3 reference type used internally
// ---------------------------------------------------------------------------

/// Minimal interface of Geodesic3 used by geodesic_line3.dart.
/// The full Geodesic3 class must implement this.
abstract class Geodesic3Ref {
  Ellipsoid3 get t;
  double get ellipthresh;
  bool get umbalt;
  GeodesicLine3? get umbline;
}

// ---------------------------------------------------------------------------
// Gamblk type (pre-computed geodesic direction parameters)
// ---------------------------------------------------------------------------

class Gamblk3 {
  final bool transpolar;
  final double gamma;
  final double gammax;
  final double kx2;
  final double kxp2;
  final double kx;
  final double kxp;
  final double nu;
  final double nup;

  const Gamblk3({
    required this.transpolar,
    required this.gamma,
    required this.gammax,
    required this.kx2,
    required this.kxp2,
    required this.kx,
    required this.kxp,
    required this.nu,
    required this.nup,
  });
}

/// Build gamblk for the equatorial (neg=false) or meridional (neg=true) limit.
Gamblk3 makeGamblkNeg(Geodesic3Ref tg, bool neg) {
  return Gamblk3(
    transpolar: neg,
    gamma: neg ? -0.0 : 0.0,
    nu: 0.0,
    nup: 1.0,
    gammax: 0.0,
    kx2: !neg ? tg.t.k2 : tg.t.kp2,
    kxp2: neg ? tg.t.k2 : tg.t.kp2,
    kx: !neg ? tg.t.k : tg.t.kp,
    kxp: neg ? tg.t.k : tg.t.kp,
  );
}

/// Build gamblk from (bet, omg, alp) Angles.
Gamblk3 makeGamblkFull(
    Geodesic3Ref tg, Angle bet, Angle omg, Angle alp) {
  final k = tg.t.k, kp = tg.t.kp;
  final k2 = tg.t.k2, kp2 = tg.t.kp2;
  final a = k * bet.cx * alp.sx, b = kp * omg.sx * alp.cx;
  var gamma = (a - b) * (a + b);
  var maxdiff = 0.0;
  if (k2 != 0.0 && kp2 != 0.0) {
    final alpdiff =
        2.0 * alp.cx * alp.sx * (k2 * _sq(bet.cx) + kp2 * _sq(omg.sx));
    final betdiff = -2.0 * bet.cx * bet.sx * k2 * _sq(alp.sx);
    final omgdiff = -2.0 * omg.cx * omg.sx * kp2 * _sq(alp.cx);
    maxdiff = math.max(alpdiff.abs(), math.max(betdiff.abs(), omgdiff.abs()));
  }
  if (gamma.abs() <= 3.0 * maxdiff * _machEps) {
    gamma = 0.0;
    if ((tg.umbalt && kp2 > 0.0) || k2 == 0.0) gamma = -gamma;
  }
  final transpolar = _signbit(gamma);
  final gammax = gamma.abs();
  final kx2f = !transpolar ? k2 : kp2;
  final kxp2f = transpolar ? k2 : kp2;
  final kxf = !transpolar ? k : kp;
  final kxpf = transpolar ? k : kp;
  final gammap = !transpolar
      ? math.sqrt(_sq(kxf * math.sqrt(_sq(bet.sx) + _sq(alp.cx * bet.cx))) +
          _sq(kxpf * omg.sx * alp.cx))
      : math.sqrt(
          _sq(kxpf * bet.cx * alp.sx) +
              _sq(kxf * math.sqrt(_sq(omg.cx) + _sq(alp.sx * omg.sx))));
  return Gamblk3(
    transpolar: transpolar,
    gamma: gamma,
    gammax: gammax,
    kx2: kx2f,
    kxp2: kxp2f,
    kx: kxf,
    kxp: kxpf,
    nu: math.sqrt(gammax) / kxf,
    nup: gammap / kxf,
  );
}

// ---------------------------------------------------------------------------
// zvals / zset — 2D bracket data structures for newt2
// ---------------------------------------------------------------------------

class ZVals3 {
  double z, fz, gz;
  ZVals3(this.z, this.fz, this.gz);
  ZVals3 copy() => ZVals3(z, fz, gz);
  bool lt(ZVals3 other) => z < other.z;
  bool eq(ZVals3 other) => z == other.z;
}

class ZSet3 {
  List<ZVals3> _s;
  ZSet3(ZVals3 lo, ZVals3 hi) : _s = [lo, hi];

  int num() => _s.length;
  ZVals3 min() => _s.first;
  ZVals3 max() => _s.last;
  ZVals3 val(int i) => _s[i];

  double bisect() {
    if (_s.length == 1) return _s[0].z;
    var maxgap = -1.0, maxind = 0;
    for (var i = 0; i < _s.length - 1; ++i) {
      final gap = _s[i + 1].z - _s[i].z;
      if (gap > maxgap) {
        maxgap = gap;
        maxind = i;
      }
    }
    return (_s[maxind].z + _s[maxind + 1].z) / 2.0;
  }

  int insert(ZVals3 t, [int flag = 0]) {
    if (t.z.isNaN) return -1;

    if (t.z < min().z) {
      t.fz = math.min(t.fz, min().fz);
      t.gz = math.min(t.gz, min().gz);
    } else if (t.z == min().z) {
      if (flag > 0) _s = [min()];
      t.z = min().z; t.fz = min().fz; t.gz = min().gz;
    } else if (t.z == max().z) {
      if (flag < 0) _s = [max()];
      t.z = max().z; t.fz = max().fz; t.gz = max().gz;
    } else if (max().z < t.z) {
      t.fz = math.max(t.fz, max().fz);
      t.gz = math.max(t.gz, max().gz);
    }
    if (!(min().z < t.z && t.z < max().z)) return -1;

    // Binary search
    var lo = 0, hi = _s.length;
    while (lo < hi) {
      final mid = (lo + hi) >> 1;
      if (_s[mid].z < t.z) lo = mid + 1;
      else hi = mid;
    }
    final p = lo;

    final ins = !(_s.length > p && _s[p].z == t.z);
    if (!ins) {
      t.z = _s[p].z; t.fz = _s[p].fz; t.gz = _s[p].gz;
    } else {
      if (p > 0 && p < _s.length) {
        t.fz = _clamp(t.fz, _s[p - 1].fz, _s[p].fz);
        t.gz = _clamp(t.gz, _s[p - 1].gz, _s[p].gz);
      }
    }

    int ind = -1;
    if (flag < 0) {
      _s = ins ? [t, ..._s.sublist(p)] : _s.sublist(p);
      if (ins) ind = 0;
    } else if (flag > 0) {
      if (ins) {
        _s = [..._s.sublist(0, p), t];
        ind = _s.length - 1;
      } else {
        _s = _s.sublist(0, p + 1);
      }
    } else if (ins) {
      _s = [..._s.sublist(0, p), t, ..._s.sublist(p)];
      ind = p;
    }
    return ind;
  }
}

void _zsetsinsert(
    ZSet3 xset, ZSet3 yset, ZVals3 xfg, ZVals3 yfg, double f0, double g0) {
  final xind = xset.insert(xfg);
  final yind = yset.insert(yfg);
  if (xind < 0 && yind < 0) return;

  var xa = [xset.min().copy(), xset.min().copy()];
  var xb = [xset.max().copy(), xset.max().copy()];
  var ya = [yset.min().copy(), yset.min().copy()];
  var yb = [yset.max().copy(), yset.max().copy()];

  for (var i = 0; i < xset.num(); ++i) {
    final xv = xset.val(i);
    for (var j = 0; j < yset.num(); ++j) {
      if (i == xind || j == yind) {
        final yv = yset.val(j);
        final f = (xv.fz - yv.fz) - f0;
        final g = (xv.gz + yv.gz) - g0;
        if (f <= 0) {
          if (g <= 0 && xa[0].z < xv.z) xa[0] = xv.copy();
          if (g >= 0 && yv.z < yb[0].z) yb[0] = yv.copy();
          if (i == xind && j == yind) {
            if (g <= 0 && xa[1].z < xv.z) xa[1] = xv.copy();
            if (g >= 0 && yv.z < yb[1].z) yb[1] = yv.copy();
          }
        }
        if (f >= 0) {
          if (g <= 0 && ya[0].z < yv.z) ya[0] = yv.copy();
          if (g >= 0 && xv.z < xb[0].z) xb[0] = xv.copy();
          if (i == xind && j == yind) {
            if (g <= 0 && ya[1].z < yv.z) ya[1] = yv.copy();
            if (g >= 0 && xv.z < xb[1].z) xb[1] = xv.copy();
          }
        }
      }
    }
  }

  var k = 0;
  for (; k < 2; ++k) {
    final f01 = (xa[k].fz - yb[k].fz) - f0;
    final f10 = (xb[k].fz - ya[k].fz) - f0;
    final g00 = (xa[k].gz + ya[k].gz) - g0;
    final g11 = (xb[k].gz + yb[k].gz) - g0;
    if (f01 <= 0 && f10 >= 0 && g00 <= 0 && g11 >= 0) break;
  }
  final kk = math.min(k, 1);
  xset.insert(xa[kk].copy(), -1);
  xset.insert(xb[kk].copy(), 1);
  yset.insert(ya[kk].copy(), -1);
  yset.insert(yb[kk].copy(), 1);
}

({double x, double y}) _zsetsbisect(ZSet3 xset, ZSet3 yset) {
  return (x: xset.bisect(), y: yset.bisect());
}

// ---------------------------------------------------------------------------
// Integrand functions
// ---------------------------------------------------------------------------

double _fthtbiax() => 1.0;
double _gthtbiax(double tht, double eps, double mu) => 0.0;

double _dfpsibiax(double s, double c, double eps, double mu) {
  final c2 = _sq(c) - mu * _sq(s);
  return eps / (1.0 + math.sqrt(1.0 - eps * c2));
}

double _gpsibiax(double s, double c, double eps, double mu) {
  final c2 = _sq(c) - mu * _sq(s);
  return math.sqrt(1.0 - eps * c2);
}

double _fthtp(double c, double kap, double kapp, double eps, double mu) {
  final c2 = kap * _sq(c);
  return math.sqrt((1.0 - eps * c2) / ((kapp + c2) * (c2 + mu)));
}

double _gthtp(double c, double kap, double kapp, double eps, double mu) {
  final c2 = kap * _sq(c);
  return c2 * math.sqrt((1.0 - eps * c2) / ((kapp + c2) * (c2 + mu)));
}

double _fup(double cn, double kap, double kapp, double eps, double mu) {
  final c2 = kap * _sq(cn);
  return math.sqrt((1.0 - eps * c2) / ((kapp + c2) * (kap + mu)));
}

double _gup(double cn, double dn, double kap, double kapp, double eps, double mu) {
  final c2 = kap * _sq(cn);
  return c2 * math.sqrt((1.0 - eps * c2) / ((kapp + c2) * (kap + mu)));
}

double _dfp(double c, double kap, double kapp, double eps) {
  final c2 = kap * _sq(c);
  final s = math.sqrt(kapp + c2);
  return eps * kap * math.sqrt(kapp) * c / (s * (1.0 + math.sqrt(1.0 - eps * c2)));
}

double _g0p(double c, double kap, double kapp, double eps) {
  final c2 = kap * _sq(c);
  return math.sqrt(kap * (1.0 - eps * c2) / (kapp + c2)) * c;
}

double _dfvp(double cn, double dn, double kap, double kapp, double eps) {
  return eps * kap * math.sqrt(kapp) * cn /
      (1.0 + math.sqrt(1.0 - eps * kap * _sq(cn)));
}

double _g0vp(double cn, double kap, double kapp, double eps) {
  final c2 = kap * _sq(cn);
  return math.sqrt(kap * (1.0 - eps * c2)) * cn;
}

double _fpsip(double s, double c, double kap, double kapp, double eps, double mu) {
  final c2 = kap * _sq(c) - mu * _sq(s);
  return math.sqrt((1.0 - eps * c2) / ((kapp + c2) * c2));
}

double _gpsip(double s, double c, double kap, double kapp, double eps, double mu) {
  final c2 = kap * _sq(c) - mu * _sq(s);
  return math.sqrt(c2 * (1.0 - eps * c2) / (kapp + c2));
}

double _fvp(double dn, double kap, double kapp, double eps, double mu) {
  final c2 = kap * _sq(dn);
  return math.sqrt((1.0 - eps * c2) / ((kapp + c2) * kap));
}

double _gvp(double cn, double dn, double kap, double kapp, double eps, double mu) {
  final dn2 = _sq(dn), c2 = kap * dn2;
  return dn2 * math.sqrt(kap * (1.0 - eps * c2) / (kapp + c2));
}

// ---------------------------------------------------------------------------
// HFun: wraps a TrigfunExt with various transformations
// ---------------------------------------------------------------------------

class HFun3 {
  final double _kap;
  final double _kapp;
  final double _eps;
  final double _mu;
  final double _sqrtmu;
  final double _sqrtkap;
  final double _sqrtkapp;
  final bool _distp;
  final bool _umb;
  final bool _meridr;
  final bool _meridl;
  final bool _biaxr;
  final bool _biaxl;
  bool _tx;
  EllipticFunction3? _ell;
  late TrigfunExt _fun;
  late double _max;

  HFun3(bool distp, double kap, double kapp, double epsParam, double mu,
      Geodesic3Ref tg)
      : _kap = kap,
        _kapp = kapp,
        _eps = epsParam,
        _mu = mu,
        _sqrtmu = math.sqrt(mu.abs()),
        _sqrtkap = math.sqrt(kap),
        _sqrtkapp = math.sqrt(kapp),
        _distp = distp,
        _umb = !tg.t.biaxial && mu == 0.0,
        _meridr = kap == 0.0 && mu == 0.0,
        _meridl = kapp == 0.0 && mu == 0.0,
        _biaxr = _biaxspecial(tg, mu) && kap == 0.0,
        _biaxl = _biaxspecial(tg, mu) && kapp == 0.0,
        _tx = false {
    final k = kap, kp = kapp, ep = epsParam;

    if (!distp) {
      if (_meridr || _biaxr) {
        _tx = false;
        _fun = TrigfunExt((_) => _fthtbiax(), _pi / 2.0, false);
      } else if (_meridl || _biaxl) {
        _tx = false;
        _fun = TrigfunExt(
            (psi) => _dfpsibiax(math.sin(psi), math.cos(psi), ep, mu),
            _pi / 2.0, false);
      } else if (mu > 0) {
        _tx = mu / (kap + mu) < tg.ellipthresh;
        if (_tx) {
          _ell = EllipticFunction3.full(
              kap / (kap + mu), 0.0, mu / (kap + mu), 1.0);
          final ellRef = _ell!;
          _fun = TrigfunExt((u) {
            final r = ellRef.sncndn(u);
            return _fup(r.cn, k, kp, ep, mu);
          }, ellRef.K, false);
        } else {
          _fun = TrigfunExt(
              (tht) => _fthtp(math.cos(tht), k, kp, ep, mu), _pi / 2.0, false);
        }
      } else if (mu < 0) {
        _tx = -mu / kap < tg.ellipthresh;
        if (_tx) {
          _ell = EllipticFunction3.full(
              (kap + mu) / kap, 0.0, -mu / kap, 1.0);
          final ellRef = _ell!;
          _fun = TrigfunExt((v) {
            final r = ellRef.sncndn(v);
            return _fvp(r.dn, k, kp, ep, mu);
          }, ellRef.K, false);
        } else {
          _fun = TrigfunExt(
              (psi) => _fpsip(math.sin(psi), math.cos(psi), k, kp, ep, mu),
              _pi / 2.0, false);
        }
      } else if (_umb) {
        _tx = kapp < tg.ellipthresh;
        if (_tx) {
          _ell = EllipticFunction3.full(kap, 0.0, kapp, 1.0);
          final ellRef = _ell!;
          _fun = TrigfunExt((v) {
            final r = ellRef.sncndn(v);
            return _dfvp(r.cn, r.dn, k, kp, ep);
          }, 2.0 * ellRef.K, true, 1.0);
        } else {
          _fun = TrigfunExt(
              (tht) => _dfp(math.cos(tht), k, kp, ep), _pi, true, 1.0);
        }
      } else {
        _tx = false;
        _fun = TrigfunExt.empty();
      }
    } else {
      // distp = true
      if (_meridr || _biaxr) {
        _tx = false;
        _fun = TrigfunExt(
            (tht) => _gthtbiax(tht, ep, mu), _pi / 2.0, false);
      } else if (_meridl || _biaxl) {
        _tx = false;
        _fun = TrigfunExt(
            (psi) => _gpsibiax(math.sin(psi), math.cos(psi), ep, mu),
            _pi / 2.0, false);
      } else if (mu > 0) {
        _tx = mu / (kap + mu) < tg.ellipthresh;
        if (_tx) {
          _ell = EllipticFunction3.full(
              kap / (kap + mu), 0.0, mu / (kap + mu), 1.0);
          final ellRef = _ell!;
          _fun = TrigfunExt((u) {
            final r = ellRef.sncndn(u);
            return _gup(r.cn, r.dn, k, kp, ep, mu);
          }, ellRef.K);
        } else {
          _fun = TrigfunExt(
              (tht) => _gthtp(math.cos(tht), k, kp, ep, mu), _pi / 2.0);
        }
      } else if (mu < 0) {
        _tx = -mu / kap < tg.ellipthresh;
        if (_tx) {
          _ell = EllipticFunction3.full(
              (kap + mu) / kap, 0.0, -mu / kap, 1.0);
          final ellRef = _ell!;
          _fun = TrigfunExt((v) {
            final r = ellRef.sncndn(v);
            return _gvp(r.cn, r.dn, k, kp, ep, mu);
          }, ellRef.K);
        } else {
          _fun = TrigfunExt(
              (psi) => _gpsip(math.sin(psi), math.cos(psi), k, kp, ep, mu),
              _pi / 2.0);
        }
      } else if (_umb) {
        _tx = kapp < tg.ellipthresh;
        if (_tx) {
          _ell = EllipticFunction3.full(kap, 0.0, kapp, 1.0);
          final ellRef = _ell!;
          _fun = TrigfunExt((v) {
            final r = ellRef.sncndn(v);
            return _g0vp(r.cn, k, kp, ep);
          }, 2.0 * ellRef.K, true);
        } else {
          _fun = TrigfunExt(
              (tht) => _g0p(math.cos(tht), k, kp, ep), _pi, true);
        }
      } else {
        _tx = false;
        _fun = TrigfunExt.empty();
      }
    }

    // Compute _max
    if (_umb) {
      _max = _fun.eval(_tx ? _ell!.K : _pi / 2.0);
    } else if (!distp) {
      _max = _biaxl
          ? _pi / 2.0 + _sqrtmu * _fun.max
          : _fun.max;
    } else {
      _max = _meridl ? _fun.eval(_pi / 2.0) : _fun.max;
    }
  }

  double eval(double u) {
    if (!_distp) {
      if (_biaxl) {
        return _modang(u, _sqrtmu) - _sqrtmu * _fun.eval(u);
      } else if (_meridl) {
        return 0.0;
      } else if (_umb) {
        final phi = _gd(u, _sqrtkapp);
        return u - _fun.eval(_tx ? _ell!.fPhi(phi) : phi);
      } else {
        return _fun.eval(u);
      }
    } else {
      if (_umb) {
        final phi = _gd(u, _sqrtkapp);
        return _fun.eval(_tx ? _ell!.fPhi(phi) : phi);
      } else {
        return _fun.eval(u);
      }
    }
  }

  double deriv(double u) {
    if (!_distp) {
      if (_biaxl) {
        return _sqrtmu / (_sq(math.cos(u)) - _mu * _sq(math.sin(u))) -
            _sqrtmu * _fun.deriv(u);
      } else if (_meridl) {
        return 0.0;
      } else if (_umb) {
        final phi = _gd(u, _sqrtkapp);
        final t = _kapp + _sq(sinh3(u));
        return 1.0 -
            _fun.deriv(_tx ? _ell!.fPhi(phi) : phi) /
                (_tx ? math.sqrt(t) : t / (_sqrtkapp * cosh3(u)));
      } else {
        return _fun.deriv(u);
      }
    } else {
      if (_umb) {
        final phi = _gd(u, _sqrtkapp);
        final t = _kapp + _sq(sinh3(u));
        return _fun.deriv(_tx ? _ell!.fPhi(phi) : phi) /
            (_tx ? math.sqrt(t) : t / (_sqrtkapp * cosh3(u)));
      } else {
        return _fun.deriv(u);
      }
    }
  }

  double fwd(double zeta) =>
      _umb ? _lam(zeta, _sqrtkapp) : (_tx ? _ell!.fPhi(zeta) : zeta);

  double rev(double w) =>
      _umb ? _gd(w, _sqrtkapp) : (_tx ? _ell!.am(w) : w);

  double get max => _max;

  double get maxPlus {
    final hp = halfPeriod;
    final sl = slope;
    return math.max(_max, hp * (sl == 0.0 ? 1.0 : sl) / 1000.0);
  }

  double get halfPeriod =>
      _umb ? double.infinity : (_tx ? _ell!.K : _pi / 2.0);

  double get slope {
    if (_umb) return 1.0;
    if (!_distp && _meridl) return 0.0;
    if (!_distp && _biaxl) return 1.0 - math.sqrt(-_mu) * _fun.slope;
    return _fun.slope;
  }

  int get nCoeffs => _fun.nCoeffs;

  double df(double u) => _fun.eval(u);
  double dfp(double u) => _fun.deriv(u);

  double root(double z, double u0,
      [List<int>? countn, List<int>? countb, double tol = 0.0]) {
    if (!_distp) {
      if (!z.isFinite) return z;
      if (_biaxl) {
        final d = max;
        final sl = slope;
        final ua = (z - d) / sl;
        final ub = (z + d) / sl;
        final u0c = math.min(ub, math.max(ua, u0));
        return Trigfun.root4(
            (u) => (eval(u), deriv(u)), z, u0c, ua, ub,
            halfPeriod, halfPeriod / sl, 1.0, countn, countb, tol);
      } else if (_umb) {
        final d2 = max.abs() + 2.0 * _machEps * math.max(1.0, z.abs());
        final ua = z - d2, ub = z + d2;
        final u0c = math.min(ub, math.max(ua, u0));
        return Trigfun.root4(
            (u) => (eval(u), deriv(u)), z, u0c, ua, ub,
            _pi / 2.0, _pi / 2.0, 1.0, countn, countb, tol);
      } else {
        return double.nan;
      }
    } else {
      if (!(z.isFinite && _umb)) return double.nan;
      if (z.abs() >= max) return _copysign(_bigValue(), z);
      final ua = -_bigValue(), ub = -ua;
      final u0c = math.min(ub, math.max(ua, u0));
      return Trigfun.root4(
          (u) => (eval(u), deriv(u)), z, u0c, ua, ub,
          _pi / 2.0, _pi / 2.0, 1.0, countn, countb, tol);
    }
  }

  double inv(double z, [List<int>? countn, List<int>? countb]) {
    if (!_distp) {
      if (_umb) return root(z, z, countn, countb);
      if (_biaxl) {
        return root(z, _modang(z / slope, 1.0 / _sqrtmu), countn, countb);
      }
      return _fun.inv1(z, countn, countb);
    } else {
      if (_biaxr) return double.nan;
      if (_umb) {
        final u0 = (math.log(_sqrtkapp / _sqrtkap *
                    math.tan(math.atan(_sqrtkap / _sqrtkapp) / max * z)) /
                (math.sqrt(1.0 - _eps * _kap) *
                    math.atan(_sqrtkap / _sqrtkapp) /
                    max));
        return root(z, u0, countn, countb);
      }
      return _fun.inv1(z, countn, countb);
    }
  }
}

// ---------------------------------------------------------------------------
// _newt2impl: 2D Newton root-finder
//   Solve: fx(x) - fy(y) = f0,  gx(x) + gy(y) = g0
// ---------------------------------------------------------------------------

({double x, double y}) _newt2Impl(
    double f0, double g0,
    HFun3 fx, HFun3 fy, HFun3 gx, HFun3 gy,
    double xa, double xb, double xscale,
    double ya, double yb, double yscale,
    double fscale, double gscale,
    [List<int>? countn, List<int>? countb]) {
  final tol = _machEps;
  final ftol = tol * fscale / 10.0;
  final gtol = tol * gscale / 10.0;
  final xtol = (math.pow(tol, 0.75) as double) * xscale;
  final ytol = (math.pow(tol, 0.75) as double) * yscale;

  var cntn = 0, cntb = 0;
  final xset = ZSet3(
      ZVals3(xa, fx.eval(xa), gx.eval(xa)),
      ZVals3(xb, fx.eval(xb), gx.eval(xb)));
  final yset = ZSet3(
      ZVals3(ya, fy.eval(ya), gy.eval(ya)),
      ZVals3(yb, fy.eval(yb), gy.eval(yb)));

  final p0 = _zsetsbisect(xset, yset);
  var x = p0.x, y = p0.y;

  var oldf = double.infinity, oldg = double.infinity;
  var olddx = double.infinity, olddy = double.infinity;
  var degen = false;
  var g0eff = g0;
  var ibis = -1;

  for (var i = 0; i < _maxit; ++i) {
    ++cntn;

    // Degenerate case: x interval is a single float and gy is constant
    if (!degen &&
        (xset.max().z - xset.min().z).abs() <=
            xset.min().z.abs() * _machEps * 2.0 &&
        yset.min().gz == yset.max().gz) {
      degen = true;
      final ga = (xset.min().gz + yset.min().gz) - g0eff;
      final gb = (xset.max().gz + yset.max().gz) - g0eff;
      if (gb < -ga) {
        x = xset.max().z;
        g0eff = xset.max().gz + yset.max().gz;
        xset.insert(xset.max().copy(), -1);
      } else {
        x = xset.min().z;
        g0eff = xset.min().gz + yset.min().gz;
        xset.insert(xset.min().copy(), 1);
      }
    }

    final xv = ZVals3(x, fx.eval(x), gx.eval(x));
    final yv = ZVals3(y, fy.eval(y), gy.eval(y));
    _zsetsinsert(xset, yset, xv, yv, f0, g0eff);

    final f = (xv.fz - yv.fz) - f0;
    final g = (xv.gz + yv.gz) - g0eff;
    if ((f.abs() <= ftol && g.abs() <= gtol) || f.isNaN || g.isNaN) break;

    final fxp = fx.deriv(x), fyp = fy.deriv(y);
    final gxp = gx.deriv(x), gyp = gy.deriv(y);
    final den = fxp * gyp + fyp * gxp;
    final dx = -(gyp * f + fyp * g) / den;
    final dy = -(-gxp * f + fxp * g) / den;
    final xn = x + dx, yn = y + dy;

    final curXa = xset.min().z, curXb = xset.max().z;
    final curYa = yset.min().z, curYb = yset.max().z;

    final cond1 = den > 0 &&
        (i < ibis + 2 ||
            ((2.0 * f.abs() < oldf || 2.0 * g.abs() < oldg) ||
                (2.0 * dx.abs() < olddx || 2.0 * dy.abs() < olddy)));
    final cond2 = xn >= curXa - xtol && xn <= curXb + xtol &&
        yn >= curYa - ytol && yn <= curYb + ytol;

    if (cond1 && cond2) {
      oldf = f.abs();
      oldg = g.abs();
      olddx = dx.abs();
      olddy = dy.abs();
      x = xn; y = yn;
      if (!(dx.abs() > xtol || dy.abs() > ytol)) break;
    } else {
      final pp = _zsetsbisect(xset, yset);
      final xnn = pp.x, ynn = pp.y;
      ++cntb;
      if (x == xnn && y == ynn) break;
      x = xnn; y = ynn;
      ibis = i;
    }
  }

  if (countn != null) countn[0] += cntn;
  if (countb != null) countb[0] += cntb;
  return (x: x, y: y);
}

({double x, double y}) _solve2(
    double f0, double g0,
    HFun3 fx, HFun3 fy, HFun3 gx, HFun3 gy,
    [List<int>? countn, List<int>? countb]) {
  final fxs = fx.slope, fys = fy.slope;
  final gxs = gx.slope, gys = gy.slope;
  final den = fxs * gys + fys * gxs;
  final qf = fx.maxPlus + fy.maxPlus;
  final qg = gx.maxPlus + gy.maxPlus;
  final Dx = (qf * gys + qg * fys) / den;
  final Dy = (qf * gxs + qg * fxs) / den;
  final x0 = (fys * g0 + gys * f0) / den;
  final y0 = (fxs * g0 - gxs * f0) / den;
  return _newt2Impl(f0, g0, fx, fy, gx, gy,
      x0 - Dx, x0 + Dx, fx.halfPeriod,
      y0 - Dy, y0 + Dy, fy.halfPeriod,
      (fx.halfPeriod * fxs + fy.halfPeriod * fys) / 2.0,
      (gx.halfPeriod * gxs + gy.halfPeriod * gys) / 2.0,
      countn, countb);
}

({double x, double y}) _solve2u(
    double d0, double s0,
    HFun3 fx, HFun3 fy, HFun3 gx, HFun3 gy,
    [List<int>? countn, List<int>? countb]) {
  final pi2 = _bigValue();
  final sbet = gx.max, somg = gy.max, stot = sbet + somg;
  final dbet = fx.max, domg = fy.max, del = dbet - domg;
  double u, v;
  if (s0.abs() - stot >= -5.0 * _machEps) {
    final t0 = _copysign(pi2, s0);
    final t1 = (d0 + (1.0 - 2.0 * (_signbit(s0) ? 1 : 0)) * del) / 2.0;
    u = t0 + t1; v = t0 - t1;
  } else if (d0.abs() > 2.0 * pi2 / 3.0 &&
      ((1.0 - 2.0 * (_signbit(d0) ? 1 : 0)) * s0 - (sbet - somg)).abs() <=
          7.0 * _machEps) {
    if (d0 > 0) { u = 2.0 * d0 / 3.0; v = -1.0 * d0 / 3.0; }
    else { u = 1.0 * d0 / 3.0; v = -2.0 * d0 / 3.0; }
  } else {
    final mm = 2.0;
    final res = _newt2Impl(d0, s0, fx, fy, gx, gy,
        -mm * pi2, mm * pi2, pi2,
        -mm * pi2, mm * pi2, pi2,
        pi2, pi2, countn, countb);
    u = res.x; v = res.y;
  }
  return (x: u, y: v);
}

// ---------------------------------------------------------------------------
// FLine3: "φ-arc" geodesic line (course function)
// ---------------------------------------------------------------------------

class FLine3 {
  final Geodesic3Ref _tg;
  final Gamblk3 _gm;
  HFun3? _fpsi;
  HFun3? _ftht;
  double _deltashift;

  FLine3(Geodesic3Ref tg, Gamblk3 gam)
      : _tg = tg,
        _gm = gam,
        _deltashift = double.nan {
    if (!gam.gamma.isNaN) {
      final sign = gam.transpolar ? -1.0 : 1.0;
      _fpsi = HFun3(false, gam.kx2, gam.kxp2,
          sign * tg.t.e2, -sign * gam.gamma, tg);
      _ftht = HFun3(false, gam.kxp2, gam.kx2,
          -sign * tg.t.e2, sign * gam.gamma, tg);
      _deltashift = gam.gamma == 0.0
          ? (tg.t.k2 > 0.0 && tg.t.kp2 > 0.0
              ? 2.0 * (_fpsi!.max - _ftht!.max)
              : 0.0)
          : double.nan;
    }
  }

  Geodesic3Ref get tg => _tg;
  Gamblk3 get gm => _gm;
  bool get transpolar => _gm.transpolar;
  double get gammax => _gm.gammax;
  double get kx => _gm.kx;
  double get kxp => _gm.kxp;
  double get kx2 => _gm.kx2;
  double get kxp2 => _gm.kxp2;
  double get nu => _gm.nu;
  double get nup => _gm.nup;
  HFun3 get fpsi => _fpsi!;
  HFun3 get ftht => _ftht!;
  double get deltashift => _deltashift;
}

// ---------------------------------------------------------------------------
// FIcs3: initial conditions for the course function
// ---------------------------------------------------------------------------

class FIcs3 {
  Angle tht1 = Angle.zero();
  Angle phi1 = Angle.zero();
  Angle alp1 = Angle.zero();
  Angle alp0 = Angle.zero();
  Angle phi0 = Angle.zero();
  Angle tht0 = Angle.zero();
  Angle psi1 = Angle.zero(); // only meaningful when f.gammax > 0
  double u0 = 0.0;
  double v0 = 0.0;
  double delta = 0.0;
  int Ex = 1;
  int Nx = 1;

  FIcs3(FLine3 f, Angle bet10, Angle omg10, Angle alp10) {
    final epsSq = _sq(_machEps);
    tht1 = omg10.copy();
    tht1.subEq(Angle.cardinal(1));
    phi1 = bet10.copy();
    alp1 = alp10.copy();
    alp0 = alp1.nearest(f.transpolar ? 2 : 1);

    if (f.transpolar) {
      final tmp = tht1; tht1 = phi1; phi1 = tmp;
      alp1.reflect(false, false, true);
    }

    if (!f.transpolar && phi1.sx == 0.0 && alp1.cx.abs() <= epsSq) {
      alp1 = Angle(alp1.sx, -epsSq, alp1.nx, true);
    }

    Ex = _signbit(alp1.sx) ? -1 : 1;
    Nx = _signbit(alp1.cx) ? -1 : 1;
    tht1 = tht1.flipsign(Ex.toDouble());

    if (f.gammax > 0.0) {
      phi0 = phi1.nearest(2);
      psi1 = Angle(
          f.kx * phi1.sx,
          phi0.cx * alp1.cx *
              math.sqrt(_sq(f.kx * phi1.cx) + _sq(f.kxp * tht1.cx)));
      u0 = f.fpsi.fwd(psi1.radians());
      v0 = f.ftht.fwd(tht1.radians());
      if (_biaxspecial(f.tg, f.gammax)) {
        delta = math.atan2(phi1.sx * alp1.sx.abs(), phi0.cx * alp1.cx) -
            math.sqrt(f.gammax) * f.fpsi.df(u0) -
            f.ftht.eval(v0);
      } else {
        delta = f.fpsi.eval(u0) - f.ftht.eval(v0);
      }
    } else if (f.gammax == 0.0) {
      if (f.kxp2 == 0.0) {
        // meridional geodesic on biaxial ellipsoid
        phi0 = phi1.nearest(2);
        tht0 = tht1.copy();
        if (phi1.cx == 0.0) {
          tht0.addEq((alp1.nearest(2).sub(alp1)).flipsign(phi1.sx * Ex));
        }
        u0 = f.fpsi.fwd(phi1.sub(phi0).radians());
        v0 = 0.0;
        delta = -f.ftht.eval(tht0.radians());
      } else {
        // umbilical direction
        if (phi1.cx.abs() < 8.0 * _machEps && tht1.cx.abs() < 8.0 * _machEps) {
          phi0 = phi1.nearest(1).add(Angle.cardinal(Nx));
          tht0 = tht1.nearest(1).add(Angle.cardinal(1));
          delta = f.deltashift / 2.0 -
              math.log((alp1.sx / alp1.cx).abs());
        } else {
          phi0 = phi1.nearest(2);
          tht0 = tht1.nearest(2);
          delta = Nx * f.fpsi.eval(_lamang(phi1.sub(phi0), f.tg.t.kp)) -
              f.ftht.eval(_lamang(tht1.sub(tht0), f.tg.t.k));
        }
      }
    }
  }

  void setquadrant(FLine3 f, int q) {
    final p = pos1(f.transpolar);
    p.alp1.setquadrant(q);
    final newFic = FIcs3(f, p.bet1, p.omg1, p.alp1);
    tht1 = newFic.tht1; phi1 = newFic.phi1;
    alp1 = newFic.alp1; psi1 = newFic.psi1;
    tht0 = newFic.tht0; phi0 = newFic.phi0;
    alp0 = newFic.alp0;
    u0 = newFic.u0; v0 = newFic.v0; delta = newFic.delta;
    Ex = newFic.Ex; Nx = newFic.Nx;
  }

  ({Angle bet1, Angle omg1, Angle alp1}) pos1(bool transpolar) {
    var bet1 = phi1.copy();
    var omg1 = tht1.flipsign(Ex.toDouble());
    var alp1 = this.alp1.copy();
    if (transpolar) {
      final tmp = bet1; bet1 = omg1; omg1 = tmp;
      alp1.reflect(false, false, true);
    }
    omg1.addEq(Angle.cardinal(1));
    return (bet1: bet1, omg1: omg1, alp1: alp1);
  }
}

// ---------------------------------------------------------------------------
// ArcPos0 on FLine3
// ---------------------------------------------------------------------------

({
  Angle bet2a,
  Angle omg2a,
  Angle alp2a,
  ({double phiw2, double thtw2, int ind2}) d
}) flineArcPos0(FLine3 f, FIcs3 fic, Angle tau12, [bool betp = true]) {
  var d = (phiw2: double.nan, thtw2: double.nan, ind2: 0);
  final psip = f.transpolar ? !betp : betp;

  var bet2a = Angle(0.0, 1.0, 0, true);
  var omg2a = Angle(0.0, 1.0, 0, true);
  var alp2a = Angle(0.0, 1.0, 0, true);

  if (f.gammax > 0.0) {
    double u2, v2;
    Angle psi2;
    if (psip) {
      psi2 = tau12.add(fic.psi1);
      v2 = f.fpsi.fwd(psi2.radians());
      u2 = f.ftht.inv(f.fpsi.eval(v2) - fic.delta);
      omg2a = Angle.fromRadians(f.ftht.rev(u2));
    } else {
      omg2a = fic.tht1.add(tau12);
      u2 = f.ftht.fwd(omg2a.radians());
      final u2x = f.ftht.eval(u2) + fic.delta;
      v2 = f.fpsi.inv(u2x);
      psi2 = Angle.fromRadians(f.fpsi.rev(v2));
    }
    bet2a = Angle(
        f.nup * psi2.sx,
        fic.phi0.cx *
            math.sqrt(_sq(psi2.cx) + _sq(f.nu * psi2.sx)),
        0, true).rebase(fic.phi0);
    final sAlp = fic.Ex *
        math.sqrt(_sq(f.kx * f.nu) + _sq(f.kxp * omg2a.cx));
    final cAlp = fic.phi0.cx * f.kx * f.nup * psi2.cx;
    alp2a = (sAlp == 0.0 && cAlp == 0.0)
        ? (f.transpolar ? Angle(1.0, 0.0) : Angle(0.0, 1.0))
        : Angle(sAlp, cAlp);
    d = (phiw2: v2, thtw2: u2, ind2: 0);

  } else if (f.gammax == 0.0) {
    double u2, v2;
    int ii;
    if (psip) {
      bet2a = fic.phi1.add(tau12.flipsign(fic.Nx.toDouble()));
      final phi2n = _remxAng(bet2a.sub(fic.phi0).flipsign(fic.Nx.toDouble()));
      u2 = f.fpsi.fwd(phi2n.val);
      final parity = (phi2n.n % 2 != 0) ? -1 : 1;
      if (f.kxp == 0.0) {
        v2 = 0.0;
        omg2a = Angle.fromRadians(-fic.Ex * fic.delta)
            .add(Angle.cardinal(2 * fic.Ex * phi2n.n));
        alp2a = fic.alp1.nearest(2)
            .add(Angle.cardinal(parity < 0 ? 2 : 0));
      } else {
        final deltax = _bigclamp(fic.delta + phi2n.n * f.deltashift, 2.0);
        v2 = f.ftht.inv(f.fpsi.eval(u2) - deltax);
        omg2a = Angle.fromRadians(parity * f.ftht.rev(v2)).rebase(fic.tht0);
        alp2a = Angle(
            f.kxp * fic.Ex * parity / _mcosh(v2, f.kx),
            fic.Nx * f.kx / _mcosh(u2, f.kxp));
        bet2a.addEq(Angle.epsilon().flipsign(fic.Nx.toDouble()));
      }
      ii = phi2n.n;
      d = (phiw2: u2, thtw2: v2, ind2: ii);
    } else {
      omg2a = fic.tht1.add(tau12);
      final tht2n = _remxAng(omg2a.sub(fic.tht0));
      v2 = f.ftht.fwd(tht2n.val);
      u2 = 0.0;
      final parity2 = (tht2n.n % 2 != 0) ? -1 : 1;
      if (f.kxp == 0.0) {
        u2 = v2 == 0.0 ? 0.0 : _copysign(_pi / 2.0, tht2n.val);
        bet2a = Angle.cardinal((v2.abs() == 0.0 ? 0 :
            _copysign(1.0, v2 * fic.Nx).round())).rebase(fic.phi0);
        alp2a = v2.abs() == 0.0
            ? Angle.cardinal(2)
            : fic.alp1.nearest(2)
                .add(Angle.cardinal(parity2 == 1 ? 0 : 2))
                .add(Angle.fromRadians(v2).flipsign(parity2 * bet2a.sx));
      } else {
        final deltax2 = _bigclamp(fic.delta + tht2n.n * f.deltashift, 2.0);
        u2 = f.fpsi.inv(f.ftht.eval(v2) + deltax2);
        final phi2 = fic.Nx * parity2 * f.fpsi.rev(u2);
        bet2a = Angle.fromRadians(phi2);
        alp2a = Angle(
            fic.Ex * f.kxp / _mcosh(v2, f.kx),
            f.kx * fic.Nx * parity2 / _mcosh(u2, f.kxp));
        omg2a.addEq(Angle.epsilon());
      }
      ii = tht2n.n;
      d = (phiw2: u2, thtw2: v2, ind2: ii);
    }
  }

  // Undo transpolar wrapping
  omg2a = omg2a.flipsign(fic.Ex.toDouble());
  if (f.transpolar) {
    final tmp = bet2a; bet2a = omg2a; omg2a = tmp;
    alp2a.reflect(false, false, true);
  }
  omg2a.addEq(Angle.cardinal(1));
  alp2a = alp2a.rebase(fic.alp0);
  return (bet2a: bet2a, omg2a: omg2a, alp2a: alp2a, d: d);
}

// ---------------------------------------------------------------------------
// Hybrid and Hybrid0 on FLine3
// ---------------------------------------------------------------------------

({
  Angle bet2a,
  Angle omg2a,
  Angle alp2a,
  ({double phiw2, double thtw2, int ind2}) d
}) flineHybrid(FLine3 f, FIcs3 fic, Angle betomg2, [bool betp = true]) {
  final psip = !f.transpolar ? betp : !betp;
  var bo2 = betomg2.copy();
  if (!betp) bo2.subEq(Angle.cardinal(1));

  late Angle tau12;
  if (psip) {
    final phi2 = bo2;
    if (f.gammax > 0.0) {
      final spsi = phi2.sx;
      final nu = f.nu, nup = f.nup;
      var cpsi = nu < nup
          ? (phi2.cx - nu) * (phi2.cx + nu)
          : (nup - phi2.sx) * (nup + phi2.sx);
      cpsi = !(cpsi > -_machEps) ? double.nan : (_signbit(cpsi) ? 0.0 : math.sqrt(cpsi));
      final angPsi = Angle(spsi, cpsi);
      var psi12 = angPsi.sub(fic.psi1).base();
      if (_signbit(psi12.sx)) {
        psi12 = Angle(0.0, _copysign(1.0, psi12.cx), 0, true);
      }
      tau12 = psi12;
    } else if (f.gammax == 0.0) {
      tau12 = (fic.Nx > 0
          ? phi2.sub(fic.phi1)
          : phi2.add(fic.phi1).add(Angle.cardinal(2))).base();
    } else {
      tau12 = Angle.nan();
    }
  } else {
    final tht2 = bo2.flipsign(fic.Ex.toDouble());
    if (f.gammax >= 0.0) {
      var tht2b = tht2.copy();
      tht2b.reflect(false, betp && fic.Ex < 0);
      var tht12 = tht2b.sub(fic.tht1);
      if (_signbit(tht12.sx)) {
        tht12 = Angle(0.0, _copysign(1.0, tht12.cx), 0, true);
      }
      tau12 = tht12;
    } else {
      tau12 = Angle.nan();
    }
  }
  return flineArcPos0(f, fic, tau12.base(), betp);
}

double flineHybrid0(FLine3 f, FIcs3 fic, Angle bet2, Angle omg2, [bool betp = true]) {
  final res = flineHybrid(f, fic, betp ? bet2 : omg2, betp);
  Ellipsoid3.angNorm(res.bet2a, res.omg2a, res.alp2a, !betp);
  if (betp) {
    res.omg2a.subEq(omg2);
    return res.omg2a.radians0();
  } else {
    res.bet2a.subEq(bet2);
    return res.bet2a.radians0();
  }
}

// ---------------------------------------------------------------------------
// GLine3: "s-arc" geodesic line (distance function)
// ---------------------------------------------------------------------------

class GLine3 {
  final Geodesic3Ref _tg;
  final Gamblk3 _gm;
  HFun3? _gpsi;
  HFun3? _gtht;
  double s0;

  GLine3(Geodesic3Ref tg, Gamblk3 gam)
      : _tg = tg,
        _gm = gam,
        s0 = 0.0 {
    if (!gam.gamma.isNaN) {
      final sign = gam.transpolar ? -1.0 : 1.0;
      _gpsi = HFun3(true, gam.kx2, gam.kxp2,
          sign * tg.t.e2, -sign * gam.gamma, tg);
      _gtht = HFun3(true, gam.kxp2, gam.kx2,
          -sign * tg.t.e2, sign * gam.gamma, tg);
      s0 = gam.gammax == 0.0 ? _gpsi!.max + _gtht!.max : 0.0;
    }
  }

  Geodesic3Ref get tg => _tg;
  Gamblk3 get gm => _gm;
  bool get transpolar => _gm.transpolar;
  double get gammax => _gm.gammax;
  double get kx => _gm.kx;
  double get kxp => _gm.kxp;
  double get kx2 => _gm.kx2;
  double get kxp2 => _gm.kxp2;
  HFun3 get gpsi => _gpsi!;
  HFun3 get gtht => _gtht!;
}

// ---------------------------------------------------------------------------
// GIcs3: initial conditions for the distance function
// ---------------------------------------------------------------------------

class GIcs3 {
  double sig1 = double.nan;
  double s13 = 0.0;

  GIcs3(GLine3 g, FIcs3 fic) {
    if (g.gammax > 0.0) {
      sig1 = g.gpsi.eval(fic.u0) + g.gtht.eval(fic.v0);
    } else if (g.gammax == 0.0) {
      if (g.kxp2 == 0.0) {
        sig1 = fic.Nx * g.gpsi.eval(fic.u0);
      } else {
        sig1 = fic.Nx * g.gpsi.eval(_lamang(fic.phi1.sub(fic.phi0), g.kxp)) +
            g.gtht.eval(_lamang(fic.tht1.sub(fic.tht0), g.kx));
      }
    }
  }
}

// Distance from GIcs3 to a point d
double glineDist(GLine3 g, GIcs3 ic,
    ({double phiw2, double thtw2, int ind2}) d) {
  final sig2 = g.gpsi.eval(d.phiw2) +
      g.gtht.eval(d.thtw2) +
      d.ind2 * 2.0 * g.s0;
  return (sig2 - ic.sig1) * g.tg.t.b;
}

// ---------------------------------------------------------------------------
// GeodesicLine3 — the public class
// ---------------------------------------------------------------------------

/// A triaxial geodesic line, produced by [Geodesic3.Line] or [Geodesic3.Inverse].
class GeodesicLine3 {
  final Geodesic3Ref _tg;
  late final FLine3 _f;
  late final FIcs3 _fic;
  late final GLine3 _g;
  late final GIcs3 _gic;

  // -------------------------------------------------------------------------
  // Constructors
  // -------------------------------------------------------------------------

  /// Umbilical-only constructor (no position — used for the umbilical cache).
  GeodesicLine3(Geodesic3Ref tg)
      : _tg = tg {
    final neg = (tg.umbalt && tg.t.kp2 > 0.0) || tg.t.k2 == 0.0;
    _f = FLine3(tg, makeGamblkNeg(tg, neg));
  }

  /// Full constructor from (bet1, omg1, alp1) Angles.
  GeodesicLine3.fromAngles(Geodesic3Ref tg, Angle bet1, Angle omg1, Angle alp1)
      : _tg = tg {
    bet1.round(); omg1.round(); alp1.round();
    final gam = makeGamblkFull(tg, bet1, omg1, alp1);
    _f = FLine3(tg, gam);
    _fic = FIcs3(_f, bet1, omg1, alp1);
    _g = GLine3(tg, gam);
    _gic = GIcs3(_g, _fic);
  }

  /// Construct from degrees.
  GeodesicLine3.fromDegrees(Geodesic3Ref tg, double lat1, double lon1, double azi1)
      : _tg = tg {
    final b1 = Angle.fromDegrees(lat1);
    final o1 = Angle.fromDegrees(lon1);
    final a1 = Angle.fromDegrees(azi1);
    b1.round(); o1.round(); a1.round();
    final gam = makeGamblkFull(tg, b1, o1, a1);
    _f = FLine3(tg, gam);
    _fic = FIcs3(_f, b1, o1, a1);
    _g = GLine3(tg, gam);
    _gic = GIcs3(_g, _fic);
  }

  /// Internal factory: build from assembled pieces (used by Geodesic3).
  GeodesicLine3.fromPieces3(
      Geodesic3Ref tg, FLine3 f, FIcs3 fic, GLine3 g, GIcs3 gic)
      : _tg = tg {
    _f = f; _fic = fic; _g = g; _gic = gic;
  }

  // -------------------------------------------------------------------------
  // Internal accessors (used by Geodesic3)
  // -------------------------------------------------------------------------

  FLine3 get fLine => _f;
  FIcs3 get fics => _fic;
  GLine3 get gLine => _g;
  GIcs3 get gics => _gic;

  // -------------------------------------------------------------------------
  // pos1 / pos1r
  // -------------------------------------------------------------------------

  ({Angle bet1, Angle omg1, Angle alp1}) pos1() =>
      _fic.pos1(_f.transpolar);

  ({double bet1, double omg1, double alp1}) pos1r() {
    final p = pos1();
    Ellipsoid3.angNorm(p.bet1, p.omg1, p.alp1);
    p.bet1.setn(); p.omg1.setn(); p.alp1.setn();
    return (
      bet1: p.bet1.degrees0(),
      omg1: p.omg1.degrees0(),
      alp1: p.alp1.degrees0()
    );
  }

  // -------------------------------------------------------------------------
  // Position(s12) — solve the direct problem
  // -------------------------------------------------------------------------

  /// Return position at arc length [s12] from the start.
  ({Angle bet2, Angle omg2, Angle alp2}) position(double s12) {
    final sig2 = _gic.sig1 + s12 / _tg.t.b;
    final f = _f, fic = _fic;
    var bet2a = Angle(0.0, 1.0, 0, true);
    var omg2a = Angle(0.0, 1.0, 0, true);
    var alp2a = Angle(0.0, 1.0, 0, true);

    if (f.gammax > 0.0) {
      final res = _solve2(fic.delta, sig2,
          f.fpsi, f.ftht, _g.gpsi, _g.gtht);
      final u2 = res.y, v2 = res.x;
      omg2a = Angle.fromRadians(f.ftht.rev(u2));
      final psi2 = Angle.fromRadians(f.fpsi.rev(v2));
      bet2a = Angle(
          f.gm.nup * psi2.sx,
          fic.phi0.cx * math.sqrt(_sq(psi2.cx) + _sq(f.gm.nu * psi2.sx)),
          0, true).rebase(fic.phi0);
      alp2a = Angle(
          fic.Ex * math.sqrt(_sq(f.kx * f.gm.nu) + _sq(f.kxp * omg2a.cx)),
          fic.phi0.cx * f.kx * f.gm.nup * psi2.cx);

    } else if (f.gammax == 0.0) {
      final sig2n = _remxReal(sig2, 2.0 * _g.s0);
      double u2, v2;

      if (f.kxp2 == 0.0) {
        // meridional
        final res2 = _solve2(fic.delta, sig2n.val,
            f.fpsi, f.ftht, _g.gpsi, _g.gtht);
        u2 = res2.x; v2 = res2.y;
        bet2a = Angle.fromRadians(u2);
        if (_signbit(bet2a.cx)) {
          bet2a = Angle(_copysign(1.0, bet2a.sx), _machEps / (1 << 11), 0, true);
        }
        final parity = (sig2n.n % 2 != 0) ? -1 : 1;
        final Ny = fic.Nx * parity;
        bet2a.reflect(_signbit(fic.phi0.cx * Ny), _signbit(fic.phi0.cx));
        bet2a = bet2a.rebase(fic.phi0);
        omg2a = Angle.fromRadians(-fic.delta)
            .add(Angle.cardinal(2 * sig2n.n));
        alp2a = Angle(fic.Ex * 0.0, fic.Nx * parity.toDouble(), 0, true);
      } else {
        var sig2nAdj = sig2n;
        if (sig2n.val - _g.s0 >= -5.0 * _machEps) {
          sig2nAdj = (val: -_g.s0, n: sig2n.n + 1);
        }
        final deltax = _bigclamp(fic.delta + sig2nAdj.n * f.deltashift, 1.0);
        final res3 = _solve2u(deltax, sig2nAdj.val,
            f.fpsi, f.ftht, _g.gpsi, _g.gtht);
        u2 = res3.x; v2 = res3.y;
        bet2a = _anglam(u2, f.kxp);
        omg2a = _anglam(v2, f.kx);
        final parity2 = (sig2nAdj.n % 2 != 0) ? -1 : 1;
        final Ny2 = fic.Nx * parity2;
        omg2a.addEq(Angle.cardinal(2 * sig2nAdj.n));
        omg2a = omg2a.add(fic.tht0);
        bet2a.reflect(_signbit(fic.phi0.cx * Ny2), _signbit(fic.phi0.cx));
        bet2a = bet2a.rebase(fic.phi0);
        alp2a = Angle(
            fic.Ex * f.kxp / _mcosh(v2, f.kx),
            f.kx * Ny2 / _mcosh(u2, f.kxp));
      }
    } else {
      bet2a = Angle.nan(); omg2a = Angle.nan(); alp2a = Angle.nan();
    }

    bet2a.round(); omg2a.round(); alp2a.round();
    omg2a = omg2a.flipsign(fic.Ex.toDouble());
    if (f.transpolar) {
      final tmp = bet2a; bet2a = omg2a; omg2a = tmp;
      alp2a.reflect(false, false, true);
    }
    alp2a = alp2a.rebase(fic.alp0);
    omg2a.addEq(Angle.cardinal(1));
    return (bet2: bet2a, omg2: omg2a, alp2: alp2a);
  }

  /// Return position at arc length [s12], result in degrees.
  ({double bet2, double omg2, double alp2}) positionDeg(double s12,
      [bool unroll = false]) {
    final r = position(s12);
    if (!unroll) {
      Ellipsoid3.angNorm(r.bet2, r.omg2, r.alp2);
      r.bet2.setn(); r.omg2.setn(); r.alp2.setn();
    }
    return (
      bet2: r.bet2.degrees0(),
      omg2: r.omg2.degrees0(),
      alp2: r.alp2.degrees0()
    );
  }
}









