/// UTM <=> Geographic (decimal degrees) and DMS <=> decimal degrees
/// conversion library for WGS84.
///
/// ==========================================================================
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
/// ### UTM conversion
/// ```dart
/// import 'package:utm_geo/utm_geo.dart';
///
/// // Decimal degrees → UTM
/// final utm = UTMConverter.fromLatLon(52.5200, 13.4050); // Berlin
/// print(utm); // UTMPoint(zone: 33N, E: 391284.107, N: 5820768.478)
///
/// // UTM → Decimal degrees
/// final geo = UTMConverter.toLatLon(utm);
/// print(geo); // GeoPoint(lat: 52.52, lon: 13.405)
/// ```
///
/// ### DMS conversion
/// ```dart
/// // Decimal → DMS
/// final dms = DMS.fromDecimal(48.8566); // (degrees: 48, minutes: 51, seconds: 23.76)
///
/// // DMS → Decimal
/// final dec = DMS.toDecimal(48, 51, 23.76, 'N'); // 48.8566...
///
/// // Format / Parse
/// print(DMS.format(48.8566));              // "48°51'23.76\"N"
/// print(DMS.parse("48°51'23.76\"N"));      // 48.8566...
/// ```
library utm_geo;

// Snyder-based UTM (~1 mm accuracy)
export 'src/utm.dart';
export 'src/dms.dart';

// Karney/Krüger 6th-order TM (~5 nm accuracy)
export 'src/krueger_tm.dart';
export 'src/karney_utm.dart';
