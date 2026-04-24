// ignore_for_file: unnecessary_non_null_assertion
// geodesic3.dart
// Dart port of Geodesic3.js / Geodesic3.{hpp,cpp}
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
import 'ellipsoid3.dart';
import 'geodesic_line3.dart';

// ignore_for_file: non_constant_identifier_names

const double _eps = 2.220446049250313e-16;
const double _pi = math.pi;
const int _maxit = 200;

double _sq(double x) => x * x;
bool _signbit(double x) =>
    x.isNegative || (x == 0.0 && (1.0 / x).isNegative);
double _bigValue() => -3.0 * math.log(_eps);

// ---------------------------------------------------------------------------
// Geodesic3 class
// ---------------------------------------------------------------------------

/// Triaxial geodesic solver.
///
/// Implements the direct and inverse geodesic problems on a triaxial ellipsoid
/// with semiaxes a ≥ b ≥ c > 0.
class Geodesic3 implements Geodesic3Ref {
  @override
  final Ellipsoid3 t;

  final bool umbalt = false;

  final double ellipthresh = math.sqrt(_eps);

  GeodesicLine3? _umbline;

  GeodesicLine3? get umbline => _umbline;

  // -------------------------------------------------------------------------
  // Constructors
  // -------------------------------------------------------------------------

  /// Construct from an [Ellipsoid3].
  Geodesic3(this.t);

  /// Construct from semiaxes [a] ≥ [b] ≥ [c] > 0.
  Geodesic3.fromAxes(double a, double b, double c)
      : t = Ellipsoid3(a, b, c);

  /// Construct from shape parameters (b, e2, k2, kp2).
  Geodesic3.fromParams(double b, double e2, double k2, double kp2)
      : t = Ellipsoid3.fromParams(b, e2, k2, kp2);

  // -------------------------------------------------------------------------
  // Inspectors (delegated to t)
  // -------------------------------------------------------------------------

  double get a => t.a;
  double get b => t.b;
  double get c => t.c;
  double get e2 => t.e2;
  double get k2 => t.k2;
  double get kp2 => t.kp2;
  double get k => t.k;
  double get kp => t.kp;
  bool get oblate => t.oblate;
  bool get prolate => t.prolate;
  bool get biaxial => t.biaxial;

  // -------------------------------------------------------------------------
  // gamma: compute gamblk from (bet, omg, alp) Angles
  // -------------------------------------------------------------------------

  Gamblk3 gamma(Angle bet, Angle omg, Angle alp) =>
      makeGamblkFull(this, bet, omg, alp);

  // -------------------------------------------------------------------------
  // Lazy-init umbilical line cache
  // -------------------------------------------------------------------------

  void _initUmbline() {
    _umbline ??= GeodesicLine3(this);
  }

  // -------------------------------------------------------------------------
  // HybridA: helper for findroot — compute longitude difference at bet2
  // -------------------------------------------------------------------------

  double hybridA(
      Angle bet1, Angle omg1, Angle alp1, Angle bet2, Angle omg2, bool betp) {
    final b1 = bet1.copy(), o1 = omg1.copy(), a1 = alp1.copy();
    final gam = makeGamblkFull(this, b1, o1, a1);
    final lf = FLine3(this, gam);
    final ic = FIcs3(lf, b1, o1, a1);
    return flineHybrid0(lf, ic, bet2, omg2, betp);
  }

  // -------------------------------------------------------------------------
  // findroot: Chandrupatla 1997 root-finding algorithm on Angle domain
  // -------------------------------------------------------------------------

  Angle findroot(
      double Function(Angle) f,
      Angle xa, Angle xb, double fa, double fb,
      [List<int>? countn, List<int>? countb]) {
    Angle? xm;
    var cntn = 0, cntb = 0;

    if (fa * fb >= 0) {
      if (fa == fb && fa.abs() <= 512 * _eps) {
        return Angle(xa.sx + xb.sx, xa.cx + xb.cx);
      } else if (math.min(fa.abs(), fb.abs()) > 2 * _eps) {
        throw StateError(
            'Geodesic3::findroot: bad inputs (same sign, large residuals)');
      } else {
        return fa.abs() < fb.abs() ? xa : xb;
      }
    }

    var trip = false;
    var t = 0.5, tp = 0.5;
    var ab = xa.sub(xb).radians0();
    var ft = 0.0, fc = fa;
    var xc = xa.copy();

    for (var iter = 0; iter < _maxit; ++iter) {
      Angle xt;
      if (2 * t == 1) {
        xt = Angle(xa.sx + xb.sx, xa.cx + xb.cx);
      } else if (t < tp) {
        xt = xa.sub(Angle.fromRadians(t * ab));
      } else {
        xt = xb.add(Angle.fromRadians(tp * ab));
      }

      if (trip) {
        xm = xt;
        break;
      }

      ++cntn;
      ft = f(xt);

      if (ft.abs() < _eps) {
        xm = xt;
        break;
      }

      if (_signbit(ft) == _signbit(fa)) {
        xc = xa.copy(); xa = xt;
        fc = fa; fa = ft;
      } else {
        xc = xb.copy(); xb = xa.copy(); xa = xt;
        fc = fb; fb = fa; fa = ft;
      }

      xm = fb.abs() < fa.abs() ? xb : xa;
      final ab2 = xa.sub(xb).radians0();
      final ca = xc.sub(xa).radians0();
      final cb = ca + ab2;
      final tl = _eps / ab2.abs();

      trip = !(2 * tl < 1);
      final tlc = math.min(1.0 / 32.0, 16.0 * tl);

      final xi = ab2 / cb;
      final xip = ca / cb;
      final phi = (fa - fb) / (fc - fb);
      final phip = (fc - fa) / (fc - fb);

      if (!trip && _sq(phip) < xip && _sq(phi) < xi) {
        t = fa / (fb - fa) * fc / (fb - fc) -
            ca / ab2 * fa / (fc - fa) * fb / (fc - fb);
        tp = fb / (fb - fa) * fc / (fc - fa) +
            cb / ab2 * fa / (fc - fa) * fb / (fc - fb);
        t = math.max(tlc, t);
        tp = math.max(tlc, tp);
      } else {
        t = tp = 0.5;
        ++cntb;
      }
      ab = ab2; // Update bracket width for next iteration's xt computation
    } // end for iter

    xm ??= fb.abs() < fa.abs() ? xb : xa;
    if (countn != null) countn[0] += cntn;
    if (countb != null) countb[0] += cntb;
    return xm;
  }

  // -------------------------------------------------------------------------
  // Inverse: solve the inverse geodesic problem
  // -------------------------------------------------------------------------

  /// Solve the inverse problem from (bet1, omg1) to (bet2, omg2).
  ///
  /// Returns the geodesic line and arc/distance information.
  ({
    GeodesicLine3 line,
    double s12,
    Angle alp1,
    Angle alp2,
    Angle bet1,
    Angle omg1,
    Angle bet2,
    Angle omg2
  }) inverse(
      dynamic bet1arg, dynamic omg1arg, dynamic bet2arg, dynamic omg2arg) {
    // Accept degrees (num) or Angle objects
    var bet1 = bet1arg is num
        ? Angle.fromDegrees(bet1arg.toDouble())
        : (bet1arg as Angle).copy();
    var omg1 = omg1arg is num
        ? Angle.fromDegrees(omg1arg.toDouble())
        : (omg1arg as Angle).copy();
    var bet2 = bet2arg is num
        ? Angle.fromDegrees(bet2arg.toDouble())
        : (bet2arg as Angle).copy();
    var omg2 = omg2arg is num
        ? Angle.fromDegrees(omg2arg.toDouble())
        : (omg2arg as Angle).copy();

    bet1.round(); omg1.round(); bet2.round(); omg2.round();

    _initUmbline();

    // Save original coordinates for final reconstruction
    final bet10 = bet1.copy(), omg10 = omg1.copy();

    // Normalize coordinates
    final flip1 = Ellipsoid3.angNorm(bet1, omg1, null, prolate);
    final flip2 = Ellipsoid3.angNorm(bet2, omg2, null, prolate);

    // Determine swap: ensure |bet1| >= |bet2| (oblate) or |omg2| >= |omg1| (prolate)
    var tmp1 = (prolate ? omg2 : bet1).copy();
    var tmp2 = (prolate ? omg1 : bet2).copy();
    tmp1.setquadrant(0); tmp2.setquadrant(0);
    var tmp12 = tmp2.sub(tmp1);
    var swap12 = tmp12.sx > 0;
    if (!biaxial && tmp12.sx == 0.0) {
      tmp1 = omg1.copy(); tmp2 = omg2.copy();
      tmp1.setquadrant(0); tmp2.setquadrant(0);
      tmp12 = tmp2.sub(tmp1);
      swap12 = tmp12.sx < 0;
    }
    if (swap12) {
      final tmpA = bet1.copy(); bet1 = bet2.copy(); bet2 = tmpA;
      final tmpB = omg1.copy(); omg1 = omg2.copy(); omg2 = tmpB;
    }

    // Rotate to canonical frame
    if (oblate) {
      omg2 = omg2.sub(omg1).base();
      omg1 = Angle.cardinal(0);
    } else if (prolate) {
      bet2 = bet2.sub(bet1.add(Angle.cardinal(1))).base();
      bet1 = Angle.cardinal(-1);
    }

    // flipz: make bet1 <= 0
    var flipz = bet1.sx > 0;
    if (flipz) {
      bet1.reflect(true);
      bet2.reflect(true);
    }

    // flipy: make omg2 >= 0 (oblate) or bet2 in [-90,90] (prolate)
    var flipy = prolate
        ? _signbit(bet2.cx)
        : (_signbit(omg1.sx) || (omg1.sx == 0.0 && _signbit(omg2.sx)));
    if (flipy) {
      if (prolate) {
        bet2.reflect(false, true);
      } else {
        omg1.reflect(true);
        omg2.reflect(true);
      }
    }

    // flipx: make omg1.c >= 0 (i.e., omg1 in [-90,90])
    var flipx = _signbit(omg1.cx);
    if (flipx) {
      omg1.reflect(false, true);
      omg2.reflect(false, true);
    }

    // flipomg: eliminate coordinate ambiguity at poles
    var flipomg = (bet2.cx == 0.0 && _signbit(omg2.sx));
    if (flipomg) omg2.reflect(true);

    var umb1 = (bet1.cx == 0.0 && omg1.sx == 0.0);
    var umb2 = (bet2.cx == 0.0 && omg2.sx == 0.0);

    final aa = k2 * _sq(bet2.cx);
    final bb = kp2 * _sq(omg2.sx);

    var lf = _umbline!.fLine;
    FIcs3? fic;
    ({double phiw2, double thtw2, int ind2}) d =
        (phiw2: double.nan, thtw2: double.nan, ind2: 0);

    double fa = double.nan, fb = double.nan;
    late Angle alpa, alpb;
    late Angle alp2Res;
    Angle? bet2a, omg2a;
    var done = false, backside = false;

    // =====================================================================
    // Case dispatch
    // =====================================================================

    if (bet1.cx * omg1.sx == 0.0 && bet2.cx * omg2.sx == 0.0) {
      // Case A.c: both points on middle ellipse
      if (umb1 && umb2 && bet2.sx > 0 && omg2.cx < 0) {
        // A.c.1: opposite umbilical points
        final alp1Umb = biaxial
            ? Angle(k, this.kp, 0, true)
            : Angle(math.exp(lf.deltashift / 2.0), 1.0);
        fic = FIcs3(lf, bet1, omg1, alp1Umb);
        final betpUmb = k2 < kp2;
        final res = flineArcPos0(lf, fic, Angle.cardinal(2), betpUmb);
        alp2Res = res.alp2a; bet2a = res.bet2a; omg2a = res.omg2a; d = res.d;
        backside = _signbit(bet2a.cx);
        done = true;
      } else if (bet1.cx == 0.0 && bet2.cx == 0.0) {
        if (bet2.sx < 0) {
          // A.c.2
          final alp1Ac2 = Angle.cardinal(oblate ? 2 : 1);
          fic = FIcs3(lf, bet1, omg1, alp1Ac2);
          final omg12Ac2 = omg2.sub(omg1);
          if (omg12Ac2.sx == 0.0 && omg12Ac2.cx < 0) {
            d = oblate
                ? (phiw2: (biaxial ? 1 : -1) * _pi / 2.0, thtw2: -_pi / 2.0, ind2: 0)
                : prolate
                    ? (phiw2: _pi / 2.0, thtw2: -_pi / 2.0, ind2: 0)
                    : (phiw2: -_bigValue(), thtw2: _bigValue(), ind2: 0);
            alp2Res = Angle.cardinal(prolate ? 1 : 0);
            bet2a = Angle(0.0, 1.0, 0, true); omg2a = Angle(0.0, 1.0, 0, true);
          } else {
            final res = flineArcPos0(lf, fic, omg12Ac2.base(), false);
            alp2Res = res.alp2a; bet2a = res.bet2a; omg2a = res.omg2a; d = res.d;
          }
          done = true;
        } else {
          // A.c.3
          final alp1Ac3 = Angle.cardinal(
              omg1.sx == 0.0 && !prolate ? 0 : -1);
          fic = FIcs3(lf, bet1, omg1, alp1Ac3);
          if (omg1.sx == 0.0 && omg2.sx == 0.0) {
            d = biaxial
                ? (phiw2: _pi / 2.0, thtw2: -_pi / 2.0, ind2: 0)
                : (phiw2: _bigValue(), thtw2: -_bigValue(), ind2: 0);
            alp2Res = Angle.cardinal(oblate ? 0 : 1);
            bet2a = Angle(0.0, 1.0, 0, true); omg2a = Angle(0.0, 1.0, 0, true);
            done = true;
          } else {
            Angle omg2aRes;
            if (omg1.sx == 0.0) {
              omg2aRes = Angle.cardinal(2);
              alp2Res = Angle(0.0, 1.0, 0, true);
              bet2a = Angle(0.0, 1.0, 0, true);
              d = (phiw2: double.nan, thtw2: double.nan, ind2: 0);
            } else {
              final res = flineArcPos0(lf, fic, Angle.cardinal(2));
              alp2Res = res.alp2a; bet2a = res.bet2a; omg2aRes = res.omg2a; d = res.d;
            }
            final diffAc3 = omg2aRes.sub(omg2);
            if (diffAc3.sx >= -_eps / 2.0) {
              final omg12Ac3 = omg2.add(omg1);
              final res = flineArcPos0(lf, fic, omg12Ac3.base(), false);
              alp2Res = res.alp2a; bet2a = res.bet2a; omg2a = res.omg2a; d = res.d;
              if (!biaxial && _signbit(omg2a.sx)) alp2Res.reflect(true);
              done = true;
            } else {
              omg2a = omg2aRes;
              alpa = Angle.cardinal(-1).add(Angle.epsilon());
              fa = diffAc3.radians0();
              final res2 = flineArcPos0(lf, fic, Angle.cardinal(-2));
              final diff2Ac3 = res2.omg2a.sub(omg2);
              alpb = alpa.neg();
              fb = diff2Ac3.radians0();
              done = false;
            }
          }
        }
      } else {
        // A.c.4
        Angle alp1Ac4;
        if (oblate) {
          alp1Ac4 = omg2.copy();
        } else if (prolate) {
          alp1Ac4 = bet2.neg();
        } else {
          int qAc4;
          if (bet1.cx == 0.0) {
            qAc4 = omg2.cx < 0 ? 1 : (omg1.sx == 0.0 && !prolate ? 0 : -1);
          } else {
            qAc4 = omg2.cx > 0 ? 0 : 2;
          }
          alp1Ac4 = Angle.cardinal(qAc4);
        }
        fic = FIcs3(lf, bet1, omg1, alp1Ac4);
        final res = prolate
            ? flineArcPos0(lf, fic, omg2.sub(omg1).base(), false)
            : flineHybrid(lf, fic, bet2.copy());
        alp2Res = res.alp2a; bet2a = res.bet2a; omg2a = res.omg2a; d = res.d;
        if (prolate && _signbit(bet2a.cx)) alp2Res.reflect(true, true);
        backside = _signbit(bet2a.cx);
        done = true;
      }
    } else if (bet1.sx == 0.0 && bet2.sx == 0.0) {
      // Case A.b: both on equator
      final omg12Ab = omg2.sub(omg1).base();
      final E = _signbit(omg12Ab.sx) ? -1 : 1;
      final alp1Ab = Angle.cardinal(E);
      final bet1Eq = bet1.copy()..reflect(true);
      lf = FLine3(this, makeGamblkFull(this, bet1Eq, omg1, alp1Ab));
      fic = FIcs3(lf, bet1Eq, omg1, alp1Ab);
      final res = flineArcPos0(lf, fic!, Angle.cardinal(2));
      bet2a = res.bet2a; omg2a = res.omg2a; alp2Res = res.alp2a; d = res.d;
      final diffAb = omg2a.sub(omg2);
      if (E * diffAb.sx >= -_eps / 2.0) {
        final res2 = flineArcPos0(lf, fic!, omg12Ab.flipsign(E.toDouble()), false);
        alp2Res = res2.alp2a; bet2a = res2.bet2a; omg2a = res2.omg2a; d = res2.d;
        done = true;
      } else {
        alpb = Angle.cardinal(-1).sub(Angle.epsilon());
        alpa = alpb.neg();
        if (E > 0) {
          fa = diffAb.radians0();
          final res2 = flineArcPos0(lf, fic!, Angle.cardinal(-2));
          final diff2Ab = res2.omg2a.sub(omg2);
          fb = diff2Ab.radians0();
        } else {
          fb = diffAb.radians0();
          final res2 = flineArcPos0(lf, fic!, Angle.cardinal(-2));
          final diff3Ab = res2.omg2a.sub(omg2);
          fa = diff3Ab.radians0();
        }
        if (fa > 0) fa -= 2 * _pi;
        if (fb < 0) fb += 2 * _pi;
        done = false;
      }
    } else if (umb1) {
      // Case B.a: umbilical point to general point
      alp2Res = Angle(kp * omg2.sx, k * bet2.cx);
      fic = FIcs3(lf, bet2, omg2, alp2Res);
      final betpBa = aa > bb;
      final deltaBa = (lf.transpolar ? -1.0 : 1.0) * fic.delta;
      final Angle alp1Ba;
      if (oblate) {
        alp1Ba = Angle.cardinal(1).sub(Angle.fromRadians(deltaBa));
      } else if (prolate) {
        alp1Ba = Angle.fromRadians(-deltaBa);
      } else {
        alp1Ba = Angle(math.exp(lf.deltashift / 2.0 - deltaBa), 1.0);
      }
      fic = FIcs3(lf, bet1, omg1, alp1Ba);
      final tau12Ba = betpBa ? bet2.sub(bet1) : omg2.sub(omg1);
      final res = flineArcPos0(lf, fic, tau12Ba.base(), betpBa);
      alp2Res = res.alp2a; bet2a = res.bet2a; omg2a = res.omg2a; d = res.d;
      done = true;
    } else if (bet1.cx == 0.0) {
      // Case B.b
      if (!_signbit(omg2.sx)) {
        alpa = Angle.cardinal(-1).add(Angle.epsilon());
        alpb = alpa.neg();
        fa = -omg2.radians();
        fb = Angle.cardinal(2).sub(omg2).radians0();
      } else {
        alpa = Angle.cardinal(1).add(Angle.epsilon());
        alpb = alpa.neg();
        fa = Angle.cardinal(2).sub(omg2).radians0();
        fb = -omg2.radians();
      }
      done = false;
    } else if (omg1.sx == 0.0) {
      // Case B.c
      if (omg2.sx > 0) {
        alpa = Angle.epsilon();
        alpb = Angle.cardinal(2).sub(alpa);
        fa = -omg2.radians0();
        fb = Angle.cardinal(2).sub(omg2).radians0();
      } else {
        alpa = Angle(-_eps / (1 << 20), -1.0, 0, true);
        alpb = Angle(-_eps / (1 << 20), 1.0, 0, true);
        fa = Angle.cardinal(2).sub(omg2).radians0();
        fb = -omg2.radians0();
      }
      done = false;
    } else {
      // Case B.d: general triaxial — bracket search
      final f4 = [0.0, 0.0, 0.0, 0.0];
      alpa = Angle(kp * omg1.sx.abs(), k * bet1.cx.abs());
      alpb = alpa.copy();
      fic = FIcs3(lf, bet1, omg1, alpb);
      var qb = 0, qa = 3;
      for (; !done && qb <= 4; ++qb, ++qa) {
        if (qb > 0) {
          alpb.setquadrant(qb);
          fic!.setquadrant(lf, qb);
        }
        if (qb < 4) {
          f4[qb & 3] = flineHybrid0(lf, fic!, bet2, omg2);
          if (f4[qb & 3].abs() < 2 * _eps) {
            final alp1Bd = alpb.copy();
            lf = FLine3(this, makeGamblkFull(this, bet1, omg1, alp1Bd));
            fic = FIcs3(lf, bet1, omg1, alp1Bd);
            final res = flineHybrid(lf, fic!, bet2.copy());
            alp2Res = res.alp2a; bet2a = res.bet2a; omg2a = res.omg2a; d = res.d;
            backside = _signbit(bet2a.cx);
            done = true;
            break;
          }
        }
        if (qb > 0 &&
            f4[qa & 3] < 0 && f4[qb & 3] > 0 &&
            f4[qb & 3] - f4[qa & 3] < 4) {
          break;
        }
      }
      if (!done) {
        fa = f4[qa & 3]; fb = f4[qb & 3];
        alpa.setquadrant(qa);
        done = false;
      }
    }

    // =====================================================================
    // Iterative search (if not done)
    // =====================================================================

    late Angle alp1;
    final countn = [0], countb = [0];

    FLine3 lfFinal;
    late FIcs3 ficFinal;
    GLine3 lg;
    GIcs3 gic;

    if (!done) {
      alp1 = findroot(
          (alp) => hybridA(bet1, omg1, alp, bet2, omg2, true),
          alpa, alpb, fa, fb, countn, countb);
      final gamF = makeGamblkFull(this, bet1, omg1, alp1);
      lfFinal = FLine3(this, gamF);
      ficFinal = FIcs3(lfFinal, bet1, omg1, alp1);
        final betpFinal = omg1.sx <= omg2.sx.abs()
          ? (lfFinal.transpolar
            ? 2.0 * (aa + lfFinal.gammax) > (aa + bb)
            : 2.0 * (bb + lfFinal.gammax) < (aa + bb))
          : true;
      final res = betpFinal
          ? flineHybrid(lfFinal, ficFinal, bet2.copy(), true)
          : flineHybrid(lfFinal, ficFinal, omg2.copy(), false);
      alp2Res = res.alp2a; bet2a = res.bet2a; omg2a = res.omg2a; d = res.d;
      backside = (betpFinal || bet2.cx != 0.0) && _signbit(bet2a.cx);
      lg = GLine3(this, gamF);
      gic = GIcs3(lg, ficFinal);
    } else {
      // done=true: extract alp1 from fic (null-safe)
      if (fic == null) {
        alp1 = Angle.cardinal(0);  // fallback: neutral azimuth
      } else {
        final p1 = fic!.pos1(lf.transpolar);
        alp1 = p1.alp1;
      }
      lfFinal = lf;
      ficFinal = fic ?? FIcs3(lf, bet1, omg1, alp1);
      lg = GLine3(this, lf.gm);
      gic = GIcs3(lg, ficFinal);
    }

    if (!biaxial && backside) alp2Res.reflect(true, true);
    alp2Res.round();

    final s12Raw = glineDist(lg, gic, d);
    gic.s13 = _signbit(s12Raw) ? 0.0 : s12Raw;

    // =====================================================================
    // Undo coordinate transforms
    // =====================================================================

    if (flipomg) {
      omg2.reflect(true);
      alp2Res.reflect(true, true);
    }
    if (flipx) {
      omg1.reflect(false, true);
      omg2.reflect(false, true);
      alp1.reflect(true);
      alp2Res.reflect(true);
    }
    if (flipy) {
      if (prolate) {
        bet2.reflect(false, true);
        alp1.reflect(false, true);
        alp2Res.reflect(false, true);
      } else {
        omg1.reflect(true);
        omg2.reflect(true);
        alp1.reflect(true);
        alp2Res.reflect(true);
      }
    }
    if (flipz) {
      bet1.reflect(true);
      bet2.reflect(true);
      alp1.reflect(false, true);
      alp2Res.reflect(false, true);
    }
    if (swap12) {
      final tmpA = bet1.copy(); bet1 = bet2.copy(); bet2 = tmpA;
      final tmpB = omg1.copy(); omg1 = omg2.copy(); omg2 = tmpB;
      final tmpC = alp1.copy(); alp1 = alp2Res.copy(); alp2Res = tmpC;
      final tmpUmb = umb1; umb1 = umb2; umb2 = tmpUmb;
      alp1.addEq(Angle.cardinal(2));
      if (umb2 && !biaxial) {
        alp2Res.addEq(Angle.cardinal(
            (_signbit(alp2Res.sx) ? -1 : 1) * bet2.sx.sign.toInt()));
      } else {
        alp2Res.addEq(Angle.cardinal(2));
      }
    }
    if (flip2) Ellipsoid3.flip(bet2, omg2, alp2Res);
    if (flip1) Ellipsoid3.flip(bet1, omg1, alp1);

    alp1.setn(); alp2Res.setn();

    // Reconstruct fic/gic from original coordinates
    final ficFinalOut = FIcs3(lfFinal, bet10, omg10, alp1);
    final gicFinalOut = GIcs3(lg, ficFinalOut);
    gicFinalOut.s13 = _signbit(s12Raw) ? 0.0 : s12Raw;

    final line = GeodesicLine3.fromPieces3(
        this, lfFinal, ficFinalOut, lg, gicFinalOut);

    return (
      line: line,
      s12: s12Raw,
      alp1: alp1,
      alp2: alp2Res,
      bet1: bet1,
      omg1: omg1,
      bet2: bet2,
      omg2: omg2
    );
  }

  /// Solve inverse problem with degree inputs/outputs.
  ({
    GeodesicLine3 line,
    double s12,
    double azi1,
    double azi2
  }) inverseDeg(double lat1, double lon1, double lat2, double lon2) {
    final r = inverse(lat1, lon1, lat2, lon2);
    r.alp1.setn(); r.alp2.setn();
    return (
      line: r.line,
      s12: r.s12,
      azi1: r.alp1.degrees0(),
      azi2: r.alp2.degrees0()
    );
  }

  // -------------------------------------------------------------------------
  // Line / Direct
  // -------------------------------------------------------------------------

  /// Create a [GeodesicLine3] from position (bet1, omg1) and heading alp1.
  GeodesicLine3 line(Angle bet1, Angle omg1, Angle alp1) =>
      GeodesicLine3.fromAngles(this, bet1, omg1, alp1);

  /// Create a [GeodesicLine3] from degrees.
  GeodesicLine3 lineDeg(double lat1, double lon1, double azi1) =>
      GeodesicLine3.fromDegrees(this, lat1, lon1, azi1);

  /// Solve the direct problem: start at (bet1, omg1) with heading alp1, travel s12.
  ({
    Angle bet2,
    Angle omg2,
    Angle alp2,
    Angle alp1,
    GeodesicLine3 line,
    double s12
  }) direct(Angle bet1, Angle omg1, Angle alp1, double s12) {
    final l = GeodesicLine3.fromAngles(this, bet1, omg1, alp1);
    final r = l.position(s12);
    return (bet2: r.bet2, omg2: r.omg2, alp2: r.alp2, alp1: alp1, line: l, s12: s12);
  }

  /// Solve the direct problem with degree inputs/outputs.
  ({
    double lat2,
    double lon2,
    double azi2,
    double azi1,
    GeodesicLine3 line,
    double s12
  }) directDeg(double lat1, double lon1, double azi1, double s12) {
    final l = GeodesicLine3.fromDegrees(this, lat1, lon1, azi1);
    final r = l.positionDeg(s12);
    return (
      lat2: r.bet2,
      lon2: r.omg2,
      azi2: r.alp2,
      azi1: azi1,
      line: l,
      s12: s12
    );
  }
}














