// utm.dart
// UTM ↔ geographic (decimal degrees) conversion — Snyder formulas
//
// Based on: Snyder (1987) "Map Projections – A Working Manual", pp. 57–64
// Ellipsoid: WGS84  (EPSG:4326 → EPSG:326xx / 327xx)
// Accuracy: ~1 mm within UTM domain
// Note: For nanometer-accuracy use KarneyUTM (karney_utm.dart) instead.
//
// ==========================================================================
// Dart Port:
//   Copyright (c) 2026 Friedhold Matz <fmatz.com@gmail.com>
//   Ported to Dart in 2026.
// Licence:
//   This Dart port is licensed under the MIT/X11 License.
// ==========================================================================

import 'dart:math';

// ---------------------------------------------------------------------------
// WGS84 constants
// ---------------------------------------------------------------------------
const double _a = 6378137.0; // semi-major axis [m]
const double _f = 1.0 / 298.257223563; // flattening
const double _k0 = 0.9996; // UTM scale factor
const double _falseEasting = 500000.0; // [m]
const double _falseNorthing = 10000000.0; // [m], southern hemisphere only

// Derived (not compile-time const in Dart because sqrt is needed)
final double _e2 = 2.0 * _f - _f * _f; // eccentricity squared
final double _ep2 = _e2 / (1.0 - _e2); // second eccentricity squared

// ---------------------------------------------------------------------------
// Data classes
// ---------------------------------------------------------------------------

/// Geographic coordinate in decimal degrees (WGS84).
class GeoPoint {
  /// Latitude in decimal degrees.  Positive = North, negative = South.
  final double latitude;

  /// Longitude in decimal degrees.  Positive = East, negative = West.
  final double longitude;

  const GeoPoint(this.latitude, this.longitude);

  @override
  String toString() => 'GeoPoint(lat: $latitude, lon: $longitude)';
}

/// UTM coordinate.
class UTMPoint {
  /// Easting in metres.
  final double easting;

  /// Northing in metres.
  final double northing;

  /// UTM zone number (1–60).
  final int zone;

  /// Hemisphere: 'N' (North) or 'S' (South).
  final String hemisphere;

  const UTMPoint({
    required this.easting,
    required this.northing,
    required this.zone,
    required this.hemisphere,
  });

  @override
  String toString() =>
      'UTMPoint(zone: $zone$hemisphere, '
      'E: ${easting.toStringAsFixed(3)}, '
      'N: ${northing.toStringAsFixed(3)})';
}

// ---------------------------------------------------------------------------
// Converter
// ---------------------------------------------------------------------------

/// Converts between geographic coordinates (decimal degrees) and UTM.
///
/// Example – forward (LatLon → UTM):
/// ```dart
/// final utm = UTMConverter.fromLatLon(52.5200, 13.4050); // Berlin
/// print(utm); // UTMPoint(zone: 33N, E: 391284.107, N: 5820768.478)
/// ```
///
/// Example – inverse (UTM → LatLon):
/// ```dart
/// final geo = UTMConverter.toLatLon(
///   UTMPoint(easting: 391284.107, northing: 5820768.478, zone: 33, hemisphere: 'N'),
/// );
/// print(geo); // GeoPoint(lat: 52.52, lon: 13.405)
/// ```
class UTMConverter {
  UTMConverter._(); // static-only class

  // -------------------------------------------------------------------------
  // Geographic → UTM
  // -------------------------------------------------------------------------

  /// Converts [lat] / [lon] (decimal degrees, WGS84) to UTM.
  ///
  /// Throws [ArgumentError] when [lat] is outside −80 … +84°.
  static UTMPoint fromLatLon(double lat, double lon) {
    if (lat < -80.0 || lat > 84.0) {
      throw ArgumentError(
          'Latitude $lat° is outside the UTM range (−80° to +84°).');
    }

    // --- Determine zone number ---
    int zone = ((lon + 180.0) / 6.0).floor() + 1;
    if (lon >= 180.0) zone = 60;

    // Special zones for Norway (zone 32) and Svalbard (zones 31/33/35/37)
    if (lat >= 56.0 && lat < 64.0 && lon >= 3.0 && lon < 12.0) zone = 32;
    if (lat >= 72.0 && lat < 84.0) {
      if (lon >= 0.0 && lon < 9.0) zone = 31;
      else if (lon >= 9.0 && lon < 21.0) zone = 33;
      else if (lon >= 21.0 && lon < 33.0) zone = 35;
      else if (lon >= 33.0 && lon < 42.0) zone = 37;
    }

    // --- Central meridian ---
    final double lambda0 = ((zone - 1) * 6 - 180 + 3) * pi / 180.0;

    // --- Convert to radians ---
    final double phi = lat * pi / 180.0;
    final double lambda = lon * pi / 180.0;

    final double sinPhi = sin(phi);
    final double cosPhi = cos(phi);
    final double tanPhi = tan(phi);

    final double N = _a / sqrt(1.0 - _e2 * sinPhi * sinPhi);
    final double T = tanPhi * tanPhi;
    final double C = _ep2 * cosPhi * cosPhi;
    final double A = cosPhi * (lambda - lambda0);

    final double M = _meridionalArc(phi);

    final double A2 = A * A;
    final double A3 = A2 * A;
    final double A4 = A3 * A;
    final double A5 = A4 * A;
    final double A6 = A5 * A;

    // Easting
    final double easting = _k0 *
            N *
            (A +
                (1.0 - T + C) * A3 / 6.0 +
                (5.0 - 18.0 * T + T * T + 72.0 * C - 58.0 * _ep2) *
                    A5 /
                    120.0) +
        _falseEasting;

    // Northing
    double northing = _k0 *
        (M +
            N *
                tanPhi *
                (A2 / 2.0 +
                    (5.0 - T + 9.0 * C + 4.0 * C * C) * A4 / 24.0 +
                    (61.0 -
                            58.0 * T +
                            T * T +
                            600.0 * C -
                            330.0 * _ep2) *
                        A6 /
                        720.0));

    final String hemisphere;
    if (lat < 0.0) {
      northing += _falseNorthing;
      hemisphere = 'S';
    } else {
      hemisphere = 'N';
    }

    return UTMPoint(
      easting: easting,
      northing: northing,
      zone: zone,
      hemisphere: hemisphere,
    );
  }

  // -------------------------------------------------------------------------
  // UTM → Geographic
  // -------------------------------------------------------------------------

  /// Converts [utm] coordinates to decimal-degree geographic coordinates.
  static GeoPoint toLatLon(UTMPoint utm) {
    final double x = utm.easting - _falseEasting;
    double y = utm.northing;
    if (utm.hemisphere == 'S') y -= _falseNorthing;

    // Central meridian
    final double lambda0 = ((utm.zone - 1) * 6 - 180 + 3) * pi / 180.0;

    // Meridional arc at this northing
    final double M = y / _k0;

    // Footpoint latitude (series from Snyder eq. 3-24 / 3-26)
    final double e2 = _e2;
    final double mu = M /
        (_a *
            (1.0 -
                e2 / 4.0 -
                3.0 * e2 * e2 / 64.0 -
                5.0 * e2 * e2 * e2 / 256.0));

    final double e1 = (1.0 - sqrt(1.0 - e2)) / (1.0 + sqrt(1.0 - e2));
    final double e12 = e1 * e1;
    final double e13 = e12 * e1;
    final double e14 = e13 * e1;

    final double phi1 = mu +
        (3.0 * e1 / 2.0 - 27.0 * e13 / 32.0) * sin(2.0 * mu) +
        (21.0 * e12 / 16.0 - 55.0 * e14 / 32.0) * sin(4.0 * mu) +
        (151.0 * e13 / 96.0) * sin(6.0 * mu) +
        (1097.0 * e14 / 512.0) * sin(8.0 * mu);

    final double sinPhi1 = sin(phi1);
    final double cosPhi1 = cos(phi1);
    final double tanPhi1 = tan(phi1);

    final double N1 = _a / sqrt(1.0 - e2 * sinPhi1 * sinPhi1);
    final double T1 = tanPhi1 * tanPhi1;
    final double C1 = _ep2 * cosPhi1 * cosPhi1;
    final double R1 =
        _a * (1.0 - e2) / pow(1.0 - e2 * sinPhi1 * sinPhi1, 1.5);

    final double D = x / (N1 * _k0);
    final double D2 = D * D;
    final double D3 = D2 * D;
    final double D4 = D3 * D;
    final double D5 = D4 * D;
    final double D6 = D5 * D;

    // Latitude
    final double lat = phi1 -
        (N1 * tanPhi1 / R1) *
            (D2 / 2.0 -
                (5.0 +
                        3.0 * T1 +
                        10.0 * C1 -
                        4.0 * C1 * C1 -
                        9.0 * _ep2) *
                    D4 /
                    24.0 +
                (61.0 +
                        90.0 * T1 +
                        298.0 * C1 +
                        45.0 * T1 * T1 -
                        252.0 * _ep2 -
                        3.0 * C1 * C1) *
                    D6 /
                    720.0);

    // Longitude
    final double lon = lambda0 +
        (D -
                (1.0 + 2.0 * T1 + C1) * D3 / 6.0 +
                (5.0 -
                        2.0 * C1 +
                        28.0 * T1 -
                        3.0 * C1 * C1 +
                        8.0 * _ep2 +
                        24.0 * T1 * T1) *
                    D5 /
                    120.0) /
            cosPhi1;

    return GeoPoint(lat * 180.0 / pi, lon * 180.0 / pi);
  }

  // -------------------------------------------------------------------------
  // Helpers
  // -------------------------------------------------------------------------

  /// Meridional arc from the equator to latitude [phi] (radians).
  /// Snyder eq. 3-21.
  static double _meridionalArc(double phi) {
    final double e2 = _e2;
    final double e4 = e2 * e2;
    final double e6 = e4 * e2;
    return _a *
        ((1.0 - e2 / 4.0 - 3.0 * e4 / 64.0 - 5.0 * e6 / 256.0) * phi -
            (3.0 * e2 / 8.0 + 3.0 * e4 / 32.0 + 45.0 * e6 / 1024.0) *
                sin(2.0 * phi) +
            (15.0 * e4 / 256.0 + 45.0 * e6 / 1024.0) * sin(4.0 * phi) -
            (35.0 * e6 / 3072.0) * sin(6.0 * phi));
  }
}
