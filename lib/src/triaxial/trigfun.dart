// trigfun.dart
// Dart port of GeographicLib::Trigfun and GeographicLib::TrigfunExt
// from Trigfun.{hpp,cpp} / Trigfun.js
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
//
// Algorithm: Fourier series representation using real FFT (radix-2 Cooley-
// Tukey).  Coefficient computation imitates kissfft's transform_real.
// The inverse Trigfun uses Newton's method + Chebfun-style standardChop.
// ==========================================================================

import 'dart:math' as math;
import 'dart:typed_data';

// ignore_for_file: non_constant_identifier_names

const double _pi = math.pi;
const double _eps = 2.220446049250313e-16;
const int _maxit = 300;

// ---------------------------------------------------------------------------
// Radix-2 in-place complex DFT (forward, no normalization)
// ---------------------------------------------------------------------------

void _fftInplace(Float64List re, Float64List im) {
  final N = re.length;
  // bit-reversal permutation
  var j = 0;
  for (var i = 1; i < N; i++) {
    var bit = N >> 1;
    for (; (j & bit) != 0; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      var t = re[i]; re[i] = re[j]; re[j] = t;
      t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }
  // Cooley-Tukey butterfly
  for (var len = 2; len <= N; len <<= 1) {
    final ang = -2.0 * _pi / len;
    final wRe0 = math.cos(ang), wIm0 = math.sin(ang);
    for (var ii = 0; ii < N; ii += len) {
      var curRe = 1.0, curIm = 0.0;
      final half = len >> 1;
      for (var k = 0; k < half; k++) {
        final p = ii + k, q = ii + k + half;
        final vRe = re[q] * curRe - im[q] * curIm;
        final vIm = re[q] * curIm + im[q] * curRe;
        re[q] = re[p] - vRe; im[q] = im[p] - vIm;
        re[p] += vRe; im[p] += vIm;
        final nRe = curRe * wRe0 - curIm * wIm0;
        curIm = curRe * wIm0 + curIm * wRe0;
        curRe = nRe;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Real FFT: H[0..M-1] → cFre[0..M/2], cFim[0..M/2]
// Matches kissfft's transform_real.
// ---------------------------------------------------------------------------

({Float64List re, Float64List im}) _transformReal(Float64List H, int M) {
  final N = M >> 1;
  final re = Float64List(N), im = Float64List(N);
  for (var j = 0; j < N; j++) {
    re[j] = H[2 * j];
    im[j] = H[2 * j + 1];
  }
  _fftInplace(re, im);

  final cFre = Float64List(N + 1), cFim = Float64List(N + 1);
  cFre[0] = re[0] + im[0]; // X[0] = DC
  cFim[0] = re[0] - im[0]; // X[N] = Nyquist (packed in imaginary)

  final halfPhi = -_pi / N;
  final twRe = math.cos(halfPhi), twIm = math.sin(halfPhi);
  var curRe = twRe, curIm = twIm;
  for (var k = 1; 2 * k < N; k++) {
    final ZkRe = re[k], ZkIm = im[k];
    final ZNkRe = re[N - k], ZNkIm = im[N - k];
    final wRe = 0.5 * (ZkRe + ZNkRe);
    final wIm = 0.5 * (ZkIm - ZNkIm);
    final zRe = 0.5 * (ZkIm + ZNkIm);
    final zIm = 0.5 * (-ZkRe + ZNkRe);
    cFre[k] = wRe + (curRe * zRe - curIm * zIm);
    cFim[k] = wIm + (curRe * zIm + curIm * zRe);
    cFre[N - k] = wRe - (curRe * zRe - curIm * zIm);
    cFim[N - k] = -(wIm - (curRe * zIm + curIm * zRe));
    final nRe = curRe * twRe - curIm * twIm;
    curIm = curRe * twIm + curIm * twRe;
    curRe = nRe;
  }
  if ((N & 1) == 0) {
    cFre[N >> 1] = re[N >> 1];
    cFim[N >> 1] = -im[N >> 1];
  }
  return (re: cFre, im: cFim);
}

// ---------------------------------------------------------------------------
// Chebfun's standardChop adapted for Trigfun
// ---------------------------------------------------------------------------

int _chop(Float64List c, [double tol = 0.0, double scale = -1.0]) {
  if (tol <= 0.0) tol = _eps;
  if (tol >= 1.0) return 1;
  final n = c.length;
  if (n < 16) return n;

  // Step 1: build normalised monotone envelope m[j]
  final m = Float64List(n);
  m[n - 1] = c[n - 1].abs();
  for (var j = n - 2; j >= 0; j--) {
    m[j] = math.max(c[j].abs(), m[j + 1]);
  }
  if (m[0] == 0.0) return 1;
  if (scale >= 0) m[0] = math.max(scale, m[0]);
  for (var j = 0; j < n; j++) {
    m[j] /= m[0];
  }

  // Step 2: find plateau point
  final logtol = math.log(tol);
  var j2 = 0, plateauPoint = n;
  for (var j = 2; j <= n; j++) {
    j2 = (1.25 * j + 5).round();
    if (j2 > n) return n;
    final e1 = m[j - 1], e2 = m[j2 - 1];
    final r = 3.0 * (1.0 - math.log(e1) / logtol);
    if (e1 == 0.0 || e2 / e1 > r) {
      plateauPoint = j - 1;
      break;
    }
  }

  // Step 3: find cutoff
  if (m[plateauPoint - 1] == 0.0) return plateauPoint;
  final tol76 = tol * math.sqrt(math.pow(tol, 1.0 / 3.0) as double);
  var j3 = 0;
  for (var j = 0; j < n; j++) {
    if (m[j] >= tol76) j3++;
  }
  if (j3 < j2) {
    j2 = j3 + 1;
    m[j2 - 1] = tol76;
  }
  final tol3 = logtol / (3.0 * (j2 - 1));
  var d = 0, best = double.infinity;
  for (var j = 0; j < j2; j++) {
    final val = math.log(m[j]) - tol3 * j;
    if (val < best) {
      best = val;
      d = j;
    }
  }
  return math.max(d, 1);
}

// ---------------------------------------------------------------------------
// Internal helper: build Trigfun from coefficients
// ---------------------------------------------------------------------------

Trigfun _trigfunFromCoeffs(
    Float64List coeff, bool odd, bool sym, double h) {
  final tf = Trigfun._empty();
  tf._coeff = Float64List.fromList(coeff);
  tf._odd = odd;
  tf._sym = sym;
  tf._h = h;
  tf._q = h / 2.0;
  tf._m = coeff.length;
  tf._n = sym ? coeff.length : coeff.length - 1;
  tf._max = -1.0;
  return tf;
}

// ---------------------------------------------------------------------------
// initbysamples: build Trigfun from function samples using FFT
// ---------------------------------------------------------------------------

Trigfun _initBysamples(
    Float64List F, bool odd, bool sym, double halfp, bool centerp) {
  final n = F.length - (!(odd || sym || centerp) ? 1 : 0);
  final M = n * (sym ? 4 : 2);
  final H = Float64List(M);
  if (!centerp) {
    if (odd) H[0] = 0.0;
    for (var i = 0; i < n; i++) {
      H[i + (odd ? 1 : 0)] = F[i];
    }
    if (!odd) H[n] = sym ? 0.0 : F[n];
    if (sym) {
      for (var i = 0; i < n; i++) {
        H[2 * n - i] = (odd ? 1.0 : -1.0) * H[i];
      }
    }
    for (var i = 1; i < M ~/ 2; i++) {
      H[M - i] = (odd ? -1.0 : 1.0) * H[i];
    }
  } else {
    for (var i = 0; i < n; i++) {
      H[i] = F[i];
    }
    if (sym) {
      for (var i = 0; i < n; i++) {
        H[2 * n - i - 1] = (odd ? 1.0 : -1.0) * H[i];
      }
    }
    for (var i = 0; i < M ~/ 2; i++) {
      H[M - i - 1] = (odd ? -1.0 : 1.0) * H[i];
    }
  }

  final cF = _transformReal(H, M);
  final cFN = cF.im[0];
  cF.im[0] = 0.0;
  cF.re[M ~/ 2] = cFN;
  cF.im[M ~/ 2] = 0.0;

  if (centerp) {
    final phStep = -_pi / M;
    for (var i = 1; i <= M ~/ 2; i++) {
      final phi = i * phStep;
      final cPhi = math.cos(phi), sPhi = math.sin(phi);
      final rr = cF.re[i], ri = cF.im[i];
      cF.re[i] = rr * cPhi - ri * sPhi;
      cF.im[i] = rr * sPhi + ri * cPhi;
    }
  }

  Float64List Hout;
  if (!sym) {
    Hout = Float64List(n + 1);
    if (!odd) {
      for (var i = 0; i <= n; i++) {
        Hout[i] = cF.re[i] / n;
      }
      Hout[0] /= 2.0;
      Hout[n] = centerp ? 0.0 : Hout[n] / 2.0;
    } else {
      for (var i = 0; i <= n; i++) {
        Hout[i] = -cF.im[i] / n;
      }
      Hout[0] = 0.0;
      Hout[n] = !centerp ? 0.0 : Hout[n] / 2.0;
    }
  } else {
    Hout = Float64List(n);
    if (!odd) {
      for (var i = 0; i < n; i++) {
        Hout[i] = cF.re[2 * i + 1] / (2 * n);
      }
    } else {
      for (var i = 0; i < n; i++) {
        Hout[i] = -cF.im[2 * i + 1] / (2 * n);
      }
    }
  }
  return _trigfunFromCoeffs(Hout, odd, sym, halfp);
}

// ---------------------------------------------------------------------------
// Trigfun class
// ---------------------------------------------------------------------------

/// Fourier approximation of a function, stored as Clenshaw coefficient array.
class Trigfun {
  bool _odd;
  bool _sym;
  double _h; // half period
  double _q; // half period / 2
  int _m; // coeff count
  int _n; // term count
  Float64List _coeff;
  double _max;

  // Default constructor: f(x) = 0
  Trigfun()
      : _odd = false,
        _sym = false,
        _h = _pi,
        _q = _pi / 2.0,
        _m = 1,
        _n = 0,
        _coeff = Float64List.fromList([0.0]),
        _max = -1.0;

  Trigfun._empty()
      : _odd = false,
        _sym = false,
        _h = _pi,
        _q = _pi / 2.0,
        _m = 1,
        _n = 0,
        _coeff = Float64List(1),
        _max = -1.0;

  // -------------------------------------------------------------------------
  // Factories
  // -------------------------------------------------------------------------

  /// Construct from [n] samples of function [f] (odd, sym, halfp, centerp).
  static Trigfun fromN(int n, double Function(double) f,
      bool odd, bool sym, double halfp,
      [bool centerp = false]) {
    if (n > 0) {
      final M = n + (!(odd || sym || centerp) ? 1 : 0);
      final p = halfp / (sym ? 2.0 : 1.0);
      final d = p / n;
      final o = centerp ? d / 2.0 : (odd ? d : 0.0);
      final F = Float64List(M);
      for (var i = 0; i < M; i++) {
        F[i] = f(o + d * i);
      }
      return _initBysamples(F, odd, sym, halfp, centerp);
    } else {
      return _trigfunFromCoeffs(
          Float64List.fromList([1.0, 0.0]), odd, sym, halfp);
    }
  }

  /// Construct from function [f] with automatic refinement.
  static Trigfun fromFunc(
      double Function(double) f, bool odd, bool sym, double halfp,
      [int nmax = 1 << 16, double tol = 0.0, double scale = -1.0]) {
    var n = 16;
    var t = fromN(n, f, odd, sym, halfp, false);
    while (n <= nmax) {
      final K = _chop(t._coeff, tol, scale < 0 ? -1.0 : scale);
      if (K < t._m) {
        t._m = K;
        t._n = t._sym ? K : K - 1;
        t._coeff = Float64List.fromList(t._coeff.sublist(0, K));
        return t;
      }
      final tx = fromN(n, f, odd, sym, halfp, true);
      t._refine(tx);
      n *= 2;
    }
    return t;
  }

  /// Construct from 2-arg function [f] with automatic refinement (for inverse).
  static Trigfun fromFunc2(
      double Function(double, double) f, bool odd, bool sym, double halfp,
      [int nmax = 1 << 16, double tol = 0.0, double scale = -1.0]) {
    var t = fromN(2, (x) => f(x, double.nan), odd, sym, halfp, false);
    while (t._n <= nmax) {
      final K = _chop(t._coeff, tol, scale < 0 ? -1.0 : scale);
      if (K < t._n) {
        t._m = K;
        t._n = t._sym ? K : K - 1;
        t._coeff = Float64List.fromList(t._coeff.sublist(0, K));
        return t;
      }
      final M = t._n;
      final p = halfp / (sym ? 2.0 : 1.0), d = p / M, o = d / 2.0;
      final F = Float64List(M);
      final tfEval = t;
      for (var i = 0; i < M; i++) {
        F[i] = f(o + d * i, tfEval.eval(o + d * i));
      }
      final tx = _initBysamples(F, odd, sym, halfp, true);
      t._refine(tx);
    }
    return t;
  }

  void _refine(Trigfun tb) {
    final m = 2 * _n + (_sym ? 0 : 1);
    final newC = Float64List(m);
    for (var i = 0; i < _n; i++) {
      newC[2 * _n + (_sym ? 0 : 1) - 1 - i] =
          (_odd ? -1.0 : 1.0) * (_coeff[i] - tb._coeff[i]) / 2.0;
    }
    if (!_sym) newC[_n] = _odd ? tb._coeff[_n] : _coeff[_n];
    for (var i = 0; i < _n; i++) {
      newC[i] = (_coeff[i] + tb._coeff[i]) / 2.0;
    }
    _max = -1.0;
    _n *= 2;
    _m = m;
    _coeff = newC;
  }

  // -------------------------------------------------------------------------
  // Evaluation (Clenshaw algorithm)
  // -------------------------------------------------------------------------

  /// Evaluate f(z) using Clenshaw summation.
  double eval(double z) {
    if (_coeff.isEmpty) return 0.0;
    final y = _pi / (_sym ? _q : _h) * z;
    var k = _m;
    final k0 = (_odd && !_sym) ? 1 : 0;
    var u0 = 0.0, u1 = 0.0;
    final x = 2.0 * math.cos(y);
    while (k > k0) {
      final tt = x * u0 - u1 + _coeff[--k];
      u1 = u0;
      u0 = tt;
    }
    if (_sym) {
      return _odd
          ? math.sin(y / 2.0) * (u0 + u1)
          : math.cos(y / 2.0) * (u0 - u1);
    } else if (!_odd) {
      return u0 - (x / 2.0) * u1;
    } else {
      return _coeff[0] * y + math.sin(y) * u0;
    }
  }

  // -------------------------------------------------------------------------
  // Properties
  // -------------------------------------------------------------------------

  /// Maximum of |coeff| (excluding constant/secular).
  double Max() {
    if (_max >= 0) return _max;
    _max = 0.0;
    for (var k = _m - 1; k >= (_sym ? 0 : 1); k--) {
      _max += _coeff[k].abs();
    }
    return _max;
  }

  double get halfPeriod => _h;
  int get nCoeffs => _m;

  double get halfRange =>
      (_odd && !_sym) ? _coeff[0] * _pi : Max();

  double get slope =>
      (_odd && !_sym) ? halfRange / _h : 0.0;

  bool get odd => _odd;
  bool get symmetric => _sym;

  void setSEcular(double f0) {
    if (!(_odd && !_sym)) {
      throw StateError('Trigfun: cannot set secular term unless odd && !sym');
    }
    _coeff[0] = f0 / _pi;
  }

  // -------------------------------------------------------------------------
  // Integral
  // -------------------------------------------------------------------------

  Trigfun integral() {
    if (_odd && !_sym && _coeff[0] != 0.0) {
      throw StateError('Trigfun: cannot integrate a secular term');
    }
    final c = Float64List.fromList(_coeff);
    final mult = (_odd ? -1.0 : 1.0) * (_sym ? _q : _h) / _pi;
    for (var i = 0; i < _m; i++) {
      c[i] *= mult / (i + (_sym ? 0.5 : 0.0));
    }
    if (!_sym) c[0] = _odd ? 0.0 : _coeff[0] * mult;
    return _trigfunFromCoeffs(c, !_odd, _sym, _h);
  }

  // -------------------------------------------------------------------------
  // Newton root finding
  // -------------------------------------------------------------------------

  /// Root finding (signature 4): general Newton with bounds.
  static double root4(
    (double, double) Function(double) ffp, // returns (f, f')
    double z,
    double x0,
    double xa,
    double xb, [
    double xscale = 1.0,
    double zscale = 1.0,
    double s = 1.0,
    List<int>? countn,
    List<int>? countb,
    double tol = 0.0,
  ]) {
    if (!(xa <= x0 && x0 <= xb)) return double.nan;
    if (x0 == xa && x0 == xb) return x0;
    final tolEff = tol <= 0 ? _eps : tol;
    final vtol = tolEff * zscale / 100.0;
    final xtol = (math.pow(tolEff, 0.75) as double) * xscale;
    var x = x0;
    var oldx = double.infinity,
        oldv = double.infinity,
        olddx = double.infinity;
    var k = 0, b = 0;
    for (; k < _maxit;) {
      k++;
      final val = ffp(x);
      final v = val.$1 - z, vp = val.$2;
      final dx = -v / vp;
      if (!(v.abs() > (k < 2 ? 0 : vtol))) break;
      if (s * v > 0) xb = math.min(xb, x);
      else xa = math.max(xa, x);
      x += dx;
      if (!(xa <= x && x <= xb) ||
          v.abs() > oldv ||
          (k > 2 && 2.0 * dx.abs() > olddx)) {
        x = (xa + xb) / 2.0;
        b++;
        if (x == oldx) break;
      } else if (!(dx.abs() > xtol)) {
        break;
      }
      oldx = x;
      oldv = v.abs();
      olddx = dx.abs();
    }
    if (countn != null) countn[0] += k;
    if (countb != null) countb[0] += b;
    return x;
  }

  /// Root finding (signature 2): for Trigfun with secular term.
  static double root2(
    Trigfun tf,
    double z,
    double Function(double) fp, [
    double x0 = double.nan,
    List<int>? countn,
    List<int>? countb,
    double tol = 0.0,
  ]) {
    if (!(tf._odd && !tf._sym)) {
      throw StateError('Trigfun: cannot take root unless odd && !sym');
    }
    final hr = _pi * tf._coeff[0];
    final sc = tf._h / hr;
    final x00 = sc * z, dx = sc.abs() * tf.Max();
    final x0eff =
        x0.isFinite ? math.min(x00 + dx, math.max(x00 - dx, x0)) : x00;
    if (dx == 0.0) return x0eff;
    final ffp = (double x) => (tf.eval(x), fp(x));
    return root4(ffp, z, x0eff, x00 - dx, x00 + dx, tf._h, hr.abs(),
        sc > 0 ? 1.0 : -1.0, countn, countb, tol);
  }

  double _inversep(
      double z, double Function(double) fp, double dx0,
      [List<int>? countn, List<int>? countb, double tol = 0.0]) {
    final hr = _pi * _coeff[0], nslope = _h / hr;
    return root2(this, z, fp, z * nslope + dx0, countn, countb, tol) -
        nslope * z;
  }

  /// Build Trigfun for f^{-1} using Newton on inversep.
  Trigfun inverse(double Function(double) fp,
      [List<int>? countn,
      List<int>? countb,
      int nmax = 1 << 16,
      double tol = 0.0,
      double scale = -1.0]) {
    if (!(_odd && !_sym && _coeff[0].isFinite && _coeff[0] != 0.0)) {
      throw StateError('Can only invert Trigfun with a secular term');
    }
    final s = _coeff[0] > 0 ? 1.0 : -1.0;
    final hp = _h, hr = _pi * _coeff[0];
    final nhp = hr * s, nhr = hp * s;
    final c0p = nhr / _pi;
    final invFunc = (double z, double dx0) =>
        _inversep(z, fp, dx0.isFinite ? dx0 : 0.0, countn, countb, tol);
    final t = fromFunc2(invFunc, _odd, _sym, nhp, nmax, tol, scale);
    t._coeff[0] = c0p;
    return t;
  }
}

// ---------------------------------------------------------------------------
// TrigfunExt class
// ---------------------------------------------------------------------------

/// Wraps a Trigfun with its derivative and optional precomputed inverse.
class TrigfunExt {
  late double Function(double) _fp;
  late bool _sym;
  late Trigfun _f;
  late double _tol;
  late int _nmax;
  bool _invp = false;
  late Trigfun _finv;

  /// Construct from derivative function [fp], half period, symmetry flag, scale.
  TrigfunExt(double Function(double) fp, double halfp,
      [bool sym = false, double scale = -1.0]) {
    _fp = fp;
    _sym = sym;
    _f = Trigfun.fromFunc(fp, false, sym, halfp, 1 << 16, 0.0, scale)
        .integral();
    _tol = math.sqrt(_eps);
    _nmax = (1.5 * _f.nCoeffs).ceil();
    _invp = false;
  }

  /// Default (empty) constructor.
  TrigfunExt.empty() {
    _fp = (_) => 0.0;
    _sym = false;
    _f = Trigfun();
    _tol = math.sqrt(_eps);
    _nmax = 1;
    _invp = false;
  }

  // -------------------------------------------------------------------------
  // Properties
  // -------------------------------------------------------------------------

  double eval(double x) => _f.eval(x);
  double deriv(double x) => _fp(x);
  int get nCoeffs => _f.nCoeffs;
  int get nCoeffsInv => _invp ? _finv.nCoeffs : -1;
  double get max => _f.Max();
  double get halfPeriod => _f.halfPeriod;
  double get halfRange => _f.halfRange;
  double get slope => _f.slope;
  bool get symmetric => _sym;

  // -------------------------------------------------------------------------
  // Inverse
  // -------------------------------------------------------------------------

  double _inv1(double z, [List<int>? countn, List<int>? countb]) {
    if (_sym) return double.nan;
    return Trigfun.root2(_f, z, _fp, double.nan, countn, countb, 0.0);
  }

  /// Public accurate inverse by direct Newton (accessible to GeodesicLine3).
  double inv1(double z, [List<int>? countn, List<int>? countb]) =>
      _inv1(z, countn, countb);

  double _inv2(double z, [List<int>? countn, List<int>? countb]) {
    if (!_invp || _sym) return double.nan;
    final x0 = _finv.eval(z);
    return Trigfun.root2(_f, z, _fp, x0, countn, countb, 0.0);
  }

  /// Compute inverse (precomputed Fourier approximation).
  void computeInverse() {
    if (!_invp && !_sym) {
      _finv = _f.inverse(_fp, null, null, _nmax, _tol, -1.0);
      _invp = true;
    }
  }

  /// Inverse f^{-1}(z): uses precomputed approximation if available.
  double inv(double z, [List<int>? countn, List<int>? countb]) {
    return _invp
        ? _inv2(z, countn, countb)
        : _inv1(z, countn, countb);
  }
}

