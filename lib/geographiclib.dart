/// Dart port of GeographicLib by Charles Karney.
///
// ==========================================================================
// Original C++: Copyright (c) Charles Karney (2008-2022)
//   <karney@alum.mit.edu>, MIT/X11 License.
//   https://geographiclib.sourceforge.io/
//   https://github.com/geographiclib/geographiclib-js
//
// Dart Port:
//   Copyright (c) 2026 Friedhold Matz <fmatz.com@gmail.com>
//   Ported to Dart in 2026.
// Licence:
//   This Dart port is licensed under the MIT/X11 License (same as original).
// ==========================================================================
/// This library re-exports all modules:
///   - Geodesic, GeodesicLine, PolygonArea  (geodesic computations)
///   - DMS                                   (degree/minute/second parsing & formatting)
///   - KruegerTM, TransverseMercatorResult   (Transverse Mercator, Krüger series)
///   - TransverseMercatorExact               (Transverse Mercator, Lee's exact method)
///   - EllipticFunction                      (Jacobi elliptic functions & integrals)
///   - KarneyUTM, UTMResult                  (UTM ↔ geographic, Karney accuracy)
///   - UTMConverter, UTMPoint, GeoPoint      (UTM ↔ geographic, Snyder formulas)
///
/// ### Quick start
/// ```dart
/// import 'package:geographiclib_dart/geographiclib.dart';
///
/// // Geodesic inverse problem
/// final g = Geodesic.WGS84;
/// final r = g.inverse(52.52, 13.41, 48.14, 11.58); // Berlin → Munich
/// print(r.s12); // distance in metres
///
/// // DMS decode/encode
/// final d = DMS.Decode('48d08\'02.7"N');
/// print(DMS.Encode(d.val, DMS.SECOND, 1, DMS.LATITUDE));
///
/// // UTM (Karney, ≤ 5 nm)
/// final utm = KarneyUTM.fromLatLon(52.52, 13.41);
/// ```
library geographiclib;

export 'src/geo_math.dart';
export 'src/geodesic.dart';
export 'src/dms.dart';
export 'src/krueger_tm.dart';
export 'src/karney_utm.dart';
export 'src/utm.dart';
export 'src/elliptic_function.dart';
export 'src/tm_exact.dart';
export 'src/exact_utm.dart';
