// karney_utm.dart
// Dart port of GeographicLib::UTMUPS (C++ → Dart)
//
// UTM ↔ geographic (decimal degrees, WGS84) converter
// using Karney's high-accuracy Krüger series (≤ 5 nm within UTM domain).
//
// Handles:
//   - UTM zone & hemisphere assignment (including Svalbard / Norway exceptions)
//   - False easting (500 000 m) and false northing (10 000 000 m south)
//   - Central meridian: lon0 = (zone − 1) × 6 − 180 + 3
//
// ==========================================================================
// Original C++: Copyright (c) Charles Karney (2008-2023)
//   <karney@alum.mit.edu>, MIT/X11 License.
//   https://geographiclib.sourceforge.io/
//   Reference: C. F. F. Karney, J. Geodesy 85:475–485 (2011)
//   https://doi.org/10.1007/s00190-011-0445-3 
//
// Dart Port:
//   Copyright (c) 2026 Friedhold Matz <fmatz.com@gmail.com>
//   Ported to Dart in 2026.
// Licence:
//   This Dart port is licensed under the MIT/X11 License (same as original).
// ==========================================================================

import 'krueger_tm.dart';

// ---------------------------------------------------------------------------
// Shared projector (WGS84, k0 = 0.9996) — constructed once
// ---------------------------------------------------------------------------
final _utm = KruegerTM(k0: 0.9996); // WGS84 defaults

// ---------------------------------------------------------------------------
// Data classes
// ---------------------------------------------------------------------------

/// UTM coordinate with full result (easting, northing, zone, hemisphere,
/// meridian convergence, and scale factor).
class KarneyUTMPoint {
  /// Easting (metres, false easting already added).
  final double easting;

  /// Northing (metres, false northing already added for southern hemisphere).
  final double northing;

  /// UTM zone number (1–60).
  final int zone;

  /// Hemisphere indicator: `'N'` or `'S'`.
  final String hemisphere;

  /// Meridian convergence (degrees), clockwise from true north.
  final double gamma;

  /// Scale factor at the point (≈ 0.9996 on central meridian).
  final double k;

  const KarneyUTMPoint({
    required this.easting,
    required this.northing,
    required this.zone,
    required this.hemisphere,
    required this.gamma,
    required this.k,
  });

  @override
  String toString() =>
      'KarneyUTMPoint(zone: $zone$hemisphere, '
      'E: ${easting.toStringAsFixed(4)}, '
      'N: ${northing.toStringAsFixed(4)}, '
      'gamma: ${gamma.toStringAsFixed(9)}°, k: ${k.toStringAsFixed(12)})';
}

/// Geographic point (decimal degrees, WGS84) with convergence and scale.
class KarneyGeoPoint {
  /// Geodetic latitude (degrees), positive north.
  final double latitude;

  /// Geodetic longitude (degrees), positive east.
  final double longitude;

  /// Meridian convergence (degrees).
  final double gamma;

  /// Scale factor.
  final double k;

  const KarneyGeoPoint({
    required this.latitude,
    required this.longitude,
    required this.gamma,
    required this.k,
  });

  @override
  String toString() =>
      'KarneyGeoPoint(lat: ${latitude.toStringAsFixed(10)}, '
      'lon: ${longitude.toStringAsFixed(10)})';
}

// ---------------------------------------------------------------------------
// Converter
// ---------------------------------------------------------------------------

/// High-accuracy UTM ↔ geographic converter based on Karney (2011).
///
/// ### Forward (geographic → UTM)
/// ```dart
/// final pt = KarneyUTMConverter.fromLatLon(52.5200, 13.4050); // Berlin
/// print(pt.easting);  // ≈ 391 779.2575 m
/// print(pt.northing); // ≈ 5 820 768.480 m
/// ```
///
/// ### Reverse (UTM → geographic)
/// ```dart
/// final geo = KarneyUTMConverter.toLatLon(pt);
/// print(geo.latitude);   // ≈ 52.520 000 000 0°
/// print(geo.longitude);  // ≈ 13.405 000 000 0°
/// ```
///
/// Accuracy: ≤ 5 nm geographic distance within the UTM domain (within 40° of
/// the central meridian).  Round-trip error is typically < 1 nm.
class KarneyUTMConverter {
  KarneyUTMConverter._();

  // ---- Zone assignment (same rules as standard UTM) -----------------------

  /// Returns the UTM zone number for [lat], [lon] (degrees).
  static int zoneFor(double lat, double lon) {
    int zone = ((lon + 180.0) / 6.0).floor() + 1;
    if (lon >= 180.0) zone = 60;

    // Norway special zone 32
    if (lat >= 56.0 && lat < 64.0 && lon >= 3.0 && lon < 12.0) zone = 32;

    // Svalbard special zones
    if (lat >= 72.0 && lat < 84.0) {
      if (lon >= 0.0 && lon < 9.0) zone = 31;
      else if (lon >= 9.0 && lon < 21.0) zone = 33;
      else if (lon >= 21.0 && lon < 33.0) zone = 35;
      else if (lon >= 33.0 && lon < 42.0) zone = 37;
    }
    return zone;
  }

  /// Returns the central meridian (degrees) for [zone].
  static double centralMeridian(int zone) => (zone - 1) * 6.0 - 180.0 + 3.0;

  // ---- Forward -----------------------------------------------------------

  /// Converts geodetic [lat] / [lon] (decimal degrees, WGS84) to UTM.
  ///
  /// Throws [ArgumentError] if [lat] is outside −80° … +84°.
  static KarneyUTMPoint fromLatLon(double lat, double lon) {
    if (lat < -80.0 || lat > 84.0) {
      throw ArgumentError(
          'Latitude $lat° is outside the UTM range (−80° to +84°).');
    }

    final int zone = zoneFor(lat, lon);
    final double lon0 = centralMeridian(zone);

    final res = _utm.forward(lon0, lat, lon);

    final bool south = lat < 0.0;
    final double easting = res.x + 500000.0;
    final double northing = res.y + (south ? 10000000.0 : 0.0);

    return KarneyUTMPoint(
      easting: easting,
      northing: northing,
      zone: zone,
      hemisphere: south ? 'S' : 'N',
      gamma: res.gamma,
      k: res.k,
    );
  }

  // ---- Reverse -----------------------------------------------------------

  /// Converts UTM [point] to geodetic coordinates (decimal degrees, WGS84).
  static KarneyGeoPoint toLatLon(KarneyUTMPoint point) {
    final double lon0 = centralMeridian(point.zone);

    final double x = point.easting - 500000.0;
    final double y = point.northing - (point.hemisphere == 'S' ? 10000000.0 : 0.0);

    final res = _utm.reverse(lon0, x, y);

    return KarneyGeoPoint(
      latitude: res.lat,
      longitude: res.lon,
      gamma: res.gamma,
      k: res.k,
    );
  }
}
