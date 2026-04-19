// tm_exact.dart
// Dart port of GeographicLib::TransverseMercatorExact (C++ → Dart)
//
// Transverse Mercator projection using Lee's (1976) exact method based on
// Jacobi elliptic functions.  Accuracy is essentially limited only by the
// precision of IEEE 754 double arithmetic (≈ 5–15 nm for WGS84).
//
// Reference:
//   L. P. Lee, "Conformal Projections Based On Jacobian Elliptic Functions",
//   Cartographica Monograph 16 (1976), Part V, pp. 67–101.
//   https://doi.org/10.3138/X687-1574-4325-WM62
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
// Notation differences from Lee (following the C++ source convention):
//   Lee x/a → xi  (northing, unit Earth)
//   Lee y/a → eta (easting,  unit Earth)
//   Lee y   → x   (easting,  metres)
//   Lee x   → y   (northing, metres)
//   Lee k   → e   (eccentricity)
//   Lee k²  → mu  (elliptic parameter)
//   Lee k'² → mv  (complementary parameter = 1 − mu)
//   Lee m   → k   (scale factor, to avoid confusion with eccentricity)

import 'dart:math' as math;
import 'elliptic_function.dart';
import 'krueger_tm.dart' show TMForwardResult, TMReverseResult;

// ---------------------------------------------------------------------------
// TransverseMercatorExact
// ---------------------------------------------------------------------------

/// Exact Transverse Mercator projection via Lee's Jacobi-elliptic method.
///
/// Provides sub-nanometre accuracy for the entire ellipsoid (not just within
/// a limited band), in contrast to the Krüger/Karney series which requires
/// the point to be within ~35° of the central meridian.
///
/// ### Usage
/// ```dart
/// // WGS84 singleton (k0 = 0.9996 for UTM)
/// final tm = TransverseMercatorExact.utm;
///
/// // Forward: geographic → projected
/// final fwd = tm.forward(lon0: 9.0, lat: 48.0, lon: 10.0);
/// print('x=${fwd.x.toStringAsFixed(3)}  y=${fwd.y.toStringAsFixed(3)}');
///
/// // Reverse: projected → geographic
/// final rev = tm.reverse(lon0: 9.0, x: fwd.x, y: fwd.y);
/// print('lat=${rev.lat}  lon=${rev.lon}');
/// ```
class TransverseMercatorExact {
  // ---- Ellipsoid parameters ------------------------------------------------
  /// Equatorial radius (metres).
  final double a;

  /// Flattening of the ellipsoid.
  final double f;

  /// Central scale factor (e.g. 0.9996 for UTM).
  final double k0;

  /// If true, extend the domain to include the entire ellipsoid.
  final bool extendp;

  // ---- Derived ---------------------------------------------------------
  late final double _mu;  // e² = f·(2−f)
  late final double _mv;  // 1 − e²
  late final double _e;   // eccentricity = √e²

  late final EllipticFunction _eEu; // parameter = μ  (k² in Lee)
  late final EllipticFunction _eEv; // parameter = mv (k'² in Lee)

  // Tolerances
  static const double _tol = 2.220446049250313e-16; // double epsilon
  static const double _tol2 = 0.1 * _tol;

  static final double _taytol = math.pow(_tol, 0.6) as double;

  // Maximum Newton iterations
  static const int _numitZeta = 6;
  static const int _numitSigma = 7;

  // ---- Constructor --------------------------------------------------------

  /// Creates an exact TM projector.
  ///
  /// [a] – equatorial radius (m); defaults to WGS84.
  /// [f] – flattening; defaults to WGS84.
  /// [k0] – central scale factor; defaults to 1.
  /// [extendp] – whether to extend the domain; defaults to false.
  TransverseMercatorExact({
    this.a = 6378137.0,
    this.f = 1.0 / 298.257223563,
    this.k0 = 1.0,
    this.extendp = false,
  }) {
    _mu = f * (2.0 - f);
    _mv = 1.0 - _mu;
    _e = math.sqrt(_mu);
    _eEu = EllipticFunction(_mu);
    _eEv = EllipticFunction(_mv);
  }

  // ---- WGS84 UTM singleton ------------------------------------------------

  /// WGS84 ellipsoid with UTM scale k0 = 0.9996.
  static final TransverseMercatorExact utm = TransverseMercatorExact(
    a: 6378137.0,
    f: 1.0 / 298.257223563,
    k0: 0.9996,
  );

  // ---- Public API ---------------------------------------------------------

  /// Forward projection: (lat, lon) in degrees → (x, y) in metres.
  ///
  /// [lon0] is the central meridian (degrees).  No false offset is applied.
  TMForwardResult forward({
    required double lon0,
    required double lat,
    required double lon,
  }) {
    lat = _latFix(lat);
    lon = _angDiff(lon0, lon);

    final int latsign = (!extendp && lat.isNegative) ? -1 : 1;
    final int lonsign = (!extendp && lon.isNegative) ? -1 : 1;
    lon *= lonsign;
    lat *= latsign;

    final bool backside = !extendp && lon > _qd;
    if (backside) {
      if (lat == 0.0) {
        // latsign unchanged (stays 1), northing will be negated via backside
      }
      lon = _hd - lon;
    }

    final double lam = lon * _degree;
    final double tau = _tand(lat);

    // u, v = Thompson TM coordinates (Lee §54)
    double u, v;
    if (lat == _qd) {
      u = _eEu.K();
      v = 0.0;
    } else if (lat == 0.0 && lon == _qd * (1.0 - _e)) {
      u = 0.0;
      v = _eEv.K();
    } else {
      _zetainvResult inv = _zetainv(_taupf(tau, _e), lam);
      u = inv.u;
      v = inv.v;
    }

    final amU = _eEu.am(u);
    final amV = _eEv.am(v);
    final snu = amU.sn; final cnu = amU.cn; final dnu = amU.dn;
    final snv = amV.sn; final cnv = amV.cn; final dnv = amV.dn;

    double xi, eta;
    final s = _sigma(snu, cnu, dnu, v, snv, cnv, dnv);
    xi = s.xi;
    eta = s.eta;

    if (backside) xi = 2.0 * _eEu.E() - xi;
    final double y = xi * a * k0 * latsign;
    final double x = eta * a * k0 * lonsign;

    double gamma, kScale;
    if (lat == _qd) {
      gamma = lon;
      kScale = 1.0;
    } else {
      // Recompute (tau, lam) from (u, v) for accurate Scale
      final z = _zeta(snu, cnu, dnu, snv, cnv, dnv);
      final tauZ = _tauf(z.taup, _e);
      final sc = _scale(tauZ, snu, cnu, dnu, snv, cnv, dnv);
      gamma = sc.gamma / _degree;
      kScale = sc.k;
    }

    if (backside) gamma = _hd - gamma;
    gamma *= latsign * lonsign;
    gamma = _angNormalize(gamma);
    kScale *= k0;

    return TMForwardResult(x: x, y: y, gamma: gamma, k: kScale);
  }

  /// Reverse projection: (x, y) in metres → (lat, lon) in degrees.
  ///
  /// [lon0] is the central meridian (degrees).  No false offset is subtracted.
  TMReverseResult reverse({
    required double lon0,
    required double x,
    required double y,
  }) {
    double xi = y / (a * k0);
    double eta = x / (a * k0);

    final int xisign = (!extendp && xi.isNegative) ? -1 : 1;
    final int etasign = (!extendp && eta.isNegative) ? -1 : 1;
    xi *= xisign;
    eta *= etasign;

    final bool backside = !extendp && xi > _eEu.E();
    if (backside) xi = 2.0 * _eEu.E() - xi;

    // u, v = Thompson TM coordinates
    double u, v;
    if (xi == 0.0 && eta == _eEv.KE()) {
      u = 0.0;
      v = _eEv.K();
    } else {
      final inv = _sigmainv(xi, eta);
      u = inv.u;
      v = inv.v;
    }

    final amU = _eEu.am(u);
    final amV = _eEv.am(v);
    final snu = amU.sn; final cnu = amU.cn; final dnu = amU.dn;
    final snv = amV.sn; final cnv = amV.cn; final dnv = amV.dn;

    double lat, lon, gamma, kScale;
    if (v != 0.0 || u != _eEu.K()) {
      final z = _zeta(snu, cnu, dnu, snv, cnv, dnv);
      final tau = _tauf(z.taup, _e);
      final phi = math.atan(tau);
      lat = phi / _degree;
      lon = z.lam / _degree;
      final sc = _scale(tau, snu, cnu, dnu, snv, cnv, dnv);
      gamma = sc.gamma / _degree;
      kScale = sc.k;
    } else {
      lat = _qd;
      lon = 0.0;
      gamma = 0.0;
      kScale = 1.0;
    }

    if (backside) lon = _hd - lon;
    lon *= etasign;
    lon = _angNormalize(lon + lon0);
    lat *= xisign;
    if (backside) gamma = _hd - gamma;
    gamma *= xisign * etasign;
    kScale *= k0;

    return TMReverseResult(lat: lat, lon: lon, gamma: gamma, k: kScale);
  }

  // ==========================================================================
  // Private algorithm routines
  // ==========================================================================

  // ---- zeta (Lee 54.17) ----------------------------------------------------

  /// Maps Thompson TM (u,v) → (taup, lam).
  ({double taup, double lam}) _zeta(
      double snu, double cnu, double dnu,
      double snv, double cnv, double dnv) {
    // Lee 54.17 rewritten to avoid overflow:
    //   atanh(snu·dnv) → asinh(snu·dnv / √(cnu² + mv·snu²·snv²))
    //   atanh(e·snu/dnv) → asinh(e·snu / √(mu·cnu² + mv·cnv²))
    const overflow = 1.0 / (2.220446049250313e-16 * 2.220446049250313e-16);
    final d1 = math.sqrt(cnu * cnu + _mv * snu * snu * snv * snv);
    final d2 = math.sqrt(_mu * cnu * cnu + _mv * cnv * cnv);
    final t1 = d1 != 0.0
        ? snu * dnv / d1
        : (snu.isNegative ? -overflow : overflow);

    // Compute t2 = sinh(e · asinh(e · snu / d2)) accurately
    final esd2 = d2 != 0.0 ? _e * snu / d2 : (snu.isNegative ? -overflow : overflow);
    final t2b = _sinh(_e * _asinh(esd2.abs())) * (esd2.isNegative ? -1.0 : 1.0);

    // psi = asinh(t1) − asinh(t2b)
    // taup = sinh(psi) = t1·√(1+t2b²) − t2b·√(1+t1²)
    final taup =
        t1 * math.sqrt(1.0 + t2b * t2b) - t2b * math.sqrt(1.0 + t1 * t1);

    final lam = (d1 != 0.0 && d2 != 0.0)
        ? math.atan2(dnu * snv, cnu * cnv) -
            _e * math.atan2(_e * cnu * snv, dnu * cnv)
        : 0.0;

    return (taup: taup, lam: lam);
  }

  // ---- dwdzeta (Lee 54.21) -------------------------------------------------

  /// Derivative dw/dzeta for Newton's method in zetainv.
  ({double du, double dv}) _dwdzeta(
      double snu, double cnu, double dnu,
      double snv, double cnv, double dnv) {
    // Lee 54.21:  (1 − dnu²·snv²) = (cnv² + mu·snu²·snv²)
    final d = _mv * _sq(_sq(cnv) + _mu * _sq(snu * snv));
    final du = cnu * dnu * dnv * (_sq(cnv) - _mu * _sq(snu * snv)) / d;
    final dv = -snu * snv * cnv * (_sq(dnu * dnv) + _mu * _sq(cnu)) / d;
    return (du: du, dv: dv);
  }

  // ---- zetainv0 (starting guess) -------------------------------------------

  /// Starting-point heuristic for [_zetainv].
  /// Returns true if the initial guess is already within tolerance.
  ({bool converged, double u, double v}) _zetainv0(double psi, double lam) {
    bool converged = false;
    double u, v;
    if (psi < -_e * math.pi / 4.0 &&
        lam > (1.0 - 2.0 * _e) * math.pi / 2.0 &&
        psi < lam - (1.0 - _e) * math.pi / 2.0) {
      // Near the south-pole log singularity
      final psix = 1.0 - psi / _e;
      final lamx = (math.pi / 2.0 - lam) / _e;
      u = _asinh(math.sin(lamx) / _hypot(math.cos(lamx), _sinh(psix))) *
          (1.0 + _mu / 2.0);
      v = math.atan2(math.cos(lamx), _sinh(psix)) * (1.0 + _mu / 2.0);
      u = _eEu.K() - u;
      v = _eEv.K() - v;
    } else if (psi < _e * math.pi / 2.0 &&
        lam > (1.0 - 2.0 * _e) * math.pi / 2.0) {
      // Near the cube-root singularity at w0 = i·Ev.K()
      final dlam = lam - (1.0 - _e) * math.pi / 2.0;
      double rad = _hypot(psi, dlam);
      final ang = math.atan2(dlam - psi, psi + dlam) - 0.75 * math.pi;
      converged = rad < _e * _taytol;
      rad = math.pow(3.0 / (_mv * _e) * rad, 1.0 / 3.0) as double;
      u = rad * math.cos(ang);
      v = rad * math.sin(ang) + _eEv.K();
    } else {
      // Spherical TM (Lee 12.6)
      v = _asinh(math.sin(lam) / _hypot(math.cos(lam), _sinh(psi)));
      u = math.atan2(_sinh(psi), math.cos(lam));
      u *= _eEu.K() / (math.pi / 2.0);
      v *= _eEu.K() / (math.pi / 2.0);
    }
    return (converged: converged, u: u, v: v);
  }

  // ---- zetainv (Newton, max 6 iterations) ---------------------------------

  _zetainvResult _zetainv(double taup, double lam) {
    final psi = _asinh(taup);
    final scal = 1.0 / _hypot(1.0, taup);
    final g0 = _zetainv0(psi, lam);
    double u = g0.u, v = g0.v;
    if (g0.converged) return _zetainvResult(u, v);
    final stol2 = _tol2 / math.max(psi * psi, 1.0);
    for (int i = 0, trip = 0; i < _numitZeta; ++i) {
      final amU = _eEu.am(u);
      final amV = _eEv.am(v);
      final snu = amU.sn; final cnu = amU.cn; final dnu = amU.dn;
      final snv = amV.sn; final cnv = amV.cn; final dnv = amV.dn;
      final z = _zeta(snu, cnu, dnu, snv, cnv, dnv);
      final dw = _dwdzeta(snu, cnu, dnu, snv, cnv, dnv);
      double tau1 = z.taup - taup;
      double lam1 = z.lam - lam;
      tau1 *= scal;
      final delu = tau1 * dw.du - lam1 * dw.dv;
      final delv = tau1 * dw.dv + lam1 * dw.du;
      u -= delu;
      v -= delv;
      if (trip != 0) break;
      final delw2 = delu * delu + delv * delv;
      if (!(delw2 >= stol2)) ++trip;
    }
    return _zetainvResult(u, v);
  }

  // ---- sigma (Lee 55.4) ---------------------------------------------------

  /// Maps Thompson TM (u,v) → GK coordinates (xi, eta).
  ({double xi, double eta}) _sigma(
      double snu, double cnu, double dnu,
      double v, double snv, double cnv, double dnv) {
    // Lee 55.4: dnu² + dnv² − 1 = mu·cnu² + mv·cnv²
    final d = _mu * _sq(cnu) + _mv * _sq(cnv);
    final xi = _eEu.eIncomplete(snu, cnu, dnu) - _mu * snu * cnu * dnu / d;
    final eta = v - _eEv.eIncomplete(snv, cnv, dnv) + _mv * snv * cnv * dnv / d;
    return (xi: xi, eta: eta);
  }

  // ---- dwdsigma (reciprocal of Lee 55.9) ----------------------------------

  /// Derivative dw/dsigma for Newton's method in sigmainv.
  ({double du, double dv}) _dwdsigma(
      double snu, double cnu, double dnu,
      double snv, double cnv, double dnv) {
    // Expanding complex dn(w)² / mv using A+S 16.21.4
    final d = _mv * _sq(_sq(cnv) + _mu * _sq(snu * snv));
    final dnr = dnu * cnv * dnv;
    final dni = -_mu * snu * cnu * snv;
    final du = (_sq(dnr) - _sq(dni)) / d;
    final dv = 2.0 * dnr * dni / d;
    return (du: du, dv: dv);
  }

  // ---- sigmainv0 (starting guess) ----------------------------------------

  ({bool converged, double u, double v}) _sigmainv0(double xi, double eta) {
    bool converged = false;
    double u, v;
    if (eta > 1.25 * _eEv.KE() ||
        (xi < -0.25 * _eEu.E() && xi < eta - _eEv.KE())) {
      // Approximate sigma as a simple pole at Eu.K() + i·Ev.K()
      final dx = xi - _eEu.E();
      final dy = eta - _eEv.KE();
      final r2 = dx * dx + dy * dy;
      u = _eEu.K() + dx / r2;
      v = _eEv.K() - dy / r2;
    } else if ((eta > 0.75 * _eEv.KE() && xi < 0.25 * _eEu.E()) ||
        eta > _eEv.KE()) {
      // Near cube-root singularity at w0 = i·Ev.K()
      final deta = eta - _eEv.KE();
      double rad = _hypot(xi, deta);
      final ang = math.atan2(deta - xi, xi + deta) - 0.75 * math.pi;
      converged = rad < 2.0 * _taytol;
      rad = math.pow(3.0 / _mv * rad, 1.0 / 3.0) as double;
      u = rad * math.cos(ang);
      v = rad * math.sin(ang) + _eEv.K();
    } else {
      // Linear approximation: w ≈ sigma · K/E
      u = xi * _eEu.K() / _eEu.E();
      v = eta * _eEu.K() / _eEu.E();
    }
    return (converged: converged, u: u, v: v);
  }

  // ---- sigmainv (Newton, max 7 iterations) --------------------------------

  _zetainvResult _sigmainv(double xi, double eta) {
    final g0 = _sigmainv0(xi, eta);
    double u = g0.u, v = g0.v;
    if (g0.converged) return _zetainvResult(u, v);
    for (int i = 0, trip = 0; i < _numitSigma; ++i) {
      final amU = _eEu.am(u);
      final amV = _eEv.am(v);
      final snu = amU.sn; final cnu = amU.cn; final dnu = amU.dn;
      final snv = amV.sn; final cnv = amV.cn; final dnv = amV.dn;
      final s = _sigma(snu, cnu, dnu, v, snv, cnv, dnv);
      final dw = _dwdsigma(snu, cnu, dnu, snv, cnv, dnv);
      final xi1 = s.xi - xi;
      final eta1 = s.eta - eta;
      final delu = xi1 * dw.du - eta1 * dw.dv;
      final delv = xi1 * dw.dv + eta1 * dw.du;
      u -= delu;
      v -= delv;
      if (trip != 0) break;
      final delw2 = delu * delu + delv * delv;
      if (!(delw2 >= _tol2 * _tol2)) ++trip;
    }
    return _zetainvResult(u, v);
  }

  // ---- Scale (Lee 55.12/55.13) --------------------------------------------

  /// Computes meridian convergence gamma (radians) and scale k at (tau, lam).
  ({double gamma, double k}) _scale(
      double tau,
      double snu, double cnu, double dnu,
      double snv, double cnv, double dnv) {
    final sec2 = 1.0 + tau * tau; // sec²(phi)
    // Lee 55.12 (negated for our sign convention)
    final gamma = math.atan2(_mv * snu * snv * cnv, cnu * dnu * dnv);
    // Lee 55.13
    final k = math.sqrt(_mv + _mu / sec2) *
        math.sqrt(sec2) *
        math.sqrt((_mv * snv * snv + cnu * cnu * dnv * dnv) /
            (_mu * cnu * cnu + _mv * cnv * cnv));
    return (gamma: gamma, k: k);
  }

  // ==========================================================================
  // Static math helpers
  // ==========================================================================

  static const double _degree = math.pi / 180.0;
  static const double _qd = 90.0;   // quarter circle in degrees
  static const double _hd = 180.0;  // half circle in degrees

  static double _sq(double x) => x * x;
  static double _hypot(double x, double y) => math.sqrt(x * x + y * y);

  static double _latFix(double x) => x.abs() > 90.0 ? double.nan : x;

  static double _angNormalize(double x) {
    x %= 360.0;
    if (x <= -180.0) return x + 360.0;
    if (x > 180.0) return x - 360.0;
    return x;
  }

  static double _angDiff(double x, double y) => _angNormalize(y - x);

  static double _tand(double x) {
    if (x.abs() == 90.0) {
      return x.isNegative ? double.negativeInfinity : double.infinity;
    }
    return math.tan(x * _degree);
  }

  static double _asinh(double x) {
    final ax = x.abs();
    if (ax > 1e7) {
      return (x < 0 ? -1.0 : 1.0) * (math.ln2 + math.log(ax));
    }
    final u = ax * ax / (math.sqrt(ax * ax + 1.0) + 1.0);
    return (x < 0 ? -1.0 : 1.0) * math.log(1.0 + ax + u);
  }

  static double _sinh(double x) {
    final e = math.exp(x);
    return (e - 1.0 / e) / 2.0;
  }

  static double _atanh(double x) {
    if (x.abs() >= 1.0) {
      return x.isNegative ? double.negativeInfinity : double.infinity;
    }
    return 0.5 * math.log((1.0 + x) / (1.0 - x));
  }

  /// τ → τ' = sinh(ψ) where ψ is the conformal latitude.
  static double _taupf(double tau, double e) {
    final tau1 = _hypot(1.0, tau);
    final sig = _sinh(e * _atanh(e * tau / tau1));
    return _hypot(1.0, sig) * tau - sig * tau1;
  }

  /// τ' → τ (Newton, 5 iterations).
  static double _tauf(double taup, double e) {
    const int numit = 5;
    final tol = math.sqrt(5e-16) / 100.0;
    final e2m = 1.0 - e * e;
    double tau = taup / e2m;
    final stol = tol * math.max(1.0, taup.abs());
    for (int i = 0; i < numit; ++i) {
      final taupa = _taupf(tau, e);
      final dtau = (taup - taupa) *
          (1.0 + e2m * tau * tau) /
          (e2m * _hypot(1.0, tau) * _hypot(1.0, taupa));
      tau += dtau;
      if (dtau.abs() < stol) break;
    }
    return tau;
  }
}

// ---------------------------------------------------------------------------
// Private result holder (avoids record allocation overhead in hot loop)
// ---------------------------------------------------------------------------
class _zetainvResult {
  final double u;
  final double v;
  _zetainvResult(this.u, this.v);
}
