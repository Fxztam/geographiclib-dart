// exact_utm.dart
// UTM ↔ geographic converter using Lee's exact Transverse Mercator.
//
// Identical zone logic and false-offset conventions as [KarneyUTMConverter],
// but backed by [TransverseMercatorExact] (Jacobi elliptic functions) instead
// of the Krüger series.  Accuracy is essentially limited only by IEEE 754
// double arithmetic (~5–15 nm) and is uniform over the full ellipsoid —
// including regions far from the central meridian or near the poles.
//
// ==========================================================================
// Original C++: Copyright (c) Charles Karney (2008-2022)
//   <karney@alum.mit.edu>, MIT/X11 License.
//   https://geographiclib.sourceforge.io/
//   Original algorithm: L. P. Lee, Cartographica Monograph 16 (1976).
//   C++ reference: GeographicLib::TransverseMercatorExact by C. Karney.
//
// Dart Port:
//   Copyright (c) 2026 Friedhold Matz <fmatz.com@gmail.com>
//   Ported to Dart in 2026.
// Licence:
//   This Dart port is licensed under the MIT/X11 License (same as original).
// ==========================================================================

import 'tm_exact.dart';
import 'karney_utm.dart' show KarneyUTMPoint, KarneyGeoPoint;

// Shared WGS84 UTM projector (constructed once).
final _tmExact = TransverseMercatorExact.utm;

// ---------------------------------------------------------------------------
// Converter
// ---------------------------------------------------------------------------

/// High-accuracy UTM ↔ geographic converter based on Lee's exact method.
///
/// Uses [TransverseMercatorExact] (Jacobi elliptic functions) instead of the
/// Krüger series, giving sub-nanometre accuracy throughout the UTM domain —
/// and beyond.  Within the normal UTM band (±3° from the central meridian),
/// results agree with [KarneyUTMConverter] to well below 1 nm.
///
/// ### Forward (geographic → UTM)
/// ```dart
/// final pt = ExactUTMConverter.fromLatLon(52.5200, 13.4050); // Berlin
/// print(pt.easting);  // ≈ 391 779.258 m
/// print(pt.northing); // ≈ 5 820 768.480 m
/// ```
///
/// ### Reverse (UTM → geographic)
/// ```dart
/// final geo = ExactUTMConverter.toLatLon(pt);
/// print(geo.latitude);   // ≈ 52.520 000 000°
/// print(geo.longitude);  // ≈ 13.405 000 000°
/// ```
///
/// Result types ([KarneyUTMPoint], [KarneyGeoPoint]) are shared with
/// [KarneyUTMConverter] so both converters are interchangeable.
class ExactUTMConverter {
  ExactUTMConverter._();

  // ---- Zone assignment (identical to KarneyUTMConverter) ------------------

  /// Returns the UTM zone number for [lat], [lon] (degrees).
  static int zoneFor(double lat, double lon) {
    int zone = ((lon + 180.0) / 6.0).floor() + 1;
    if (lon >= 180.0) zone = 60;

    // Norway special zone 32
    if (lat >= 56.0 && lat < 64.0 && lon >= 3.0 && lon < 12.0) zone = 32;

    // Svalbard special zones
    if (lat >= 72.0 && lat < 84.0) {
      if      (lon >=  0.0 && lon <  9.0) zone = 31;
      else if (lon >=  9.0 && lon < 21.0) zone = 33;
      else if (lon >= 21.0 && lon < 33.0) zone = 35;
      else if (lon >= 33.0 && lon < 42.0) zone = 37;
    }
    return zone;
  }

  /// Returns the central meridian (degrees) for [zone].
  static double centralMeridian(int zone) => (zone - 1) * 6.0 - 180.0 + 3.0;

  // ---- Forward ------------------------------------------------------------

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

    final res = _tmExact.forward(lon0: lon0, lat: lat, lon: lon);

    final bool south = lat < 0.0;
    final double easting  = res.x + 500000.0;
    final double northing = res.y + (south ? 10000000.0 : 0.0);

    return KarneyUTMPoint(
      easting:    easting,
      northing:   northing,
      zone:       zone,
      hemisphere: south ? 'S' : 'N',
      gamma:      res.gamma,
      k:          res.k,
    );
  }

  // ---- Reverse ------------------------------------------------------------

  /// Converts UTM [point] to geodetic coordinates (decimal degrees, WGS84).
  static KarneyGeoPoint toLatLon(KarneyUTMPoint point) {
    final double lon0 = centralMeridian(point.zone);

    final double x = point.easting  - 500000.0;
    final double y = point.northing - (point.hemisphere == 'S' ? 10000000.0 : 0.0);

    final res = _tmExact.reverse(lon0: lon0, x: x, y: y);

    return KarneyGeoPoint(
      latitude:  res.lat,
      longitude: res.lon,
      gamma:     res.gamma,
      k:         res.k,
    );
  }
}
