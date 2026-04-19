// test-geolib.dart
// Demo program for the geographiclib_dart package.
//
// Run:
//   dart pub get
//   dart run bin/test-geolib.dart
//
// Expected output:
//   === Geodesic inverse: Berlin → Munich ===
//   Distance s12 = 504289.453 m
//   Azimuth  azi1 = -164.317446°
//
//   === DMS decode / encode ===
//   Decoded value  = 48.134083333333336°
//   Re-encoded DMS = 48°08'02.7"N
//
//   === UTM Karney (Berlin 52.52°N 13.41°E) ===
//   Zone:      33N
//   Easting:   392118.487 m
//   Northing:  5820064.675 m
//   γ (conv):  -1.261890431°
//   k (scale): 0.999742858842
//
//   === KruegerTM forward (lon0=9°, lat=52°, lon=13.41°) ===
//   x (easting from CM):  302797.714552 m
//   y (northing):         5772534.263476 m
//   γ (convergence):      3.477748932°
//   k (scale):            1.001125279691

import 'package:geographiclib_dart/geographiclib.dart';

void main() {
  // ── 1. Geodesic inverse: distance Berlin → Munich ─────────────────────────
  final g = Geodesic.WGS84;
  final r = g.inverse(52.52, 13.41, 48.14, 11.58);
  print('=== Geodesic inverse: Berlin → Munich ===');
  print('Distance s12 = ${r.s12!.toStringAsFixed(3)} m');   // 504289.453 m
  print('Azimuth  azi1 = ${r.azi1!.toStringAsFixed(6)}°');  // -164.317446°
  print('');

  // ── 2. DMS decode / encode ─────────────────────────────────────────────────
  print('=== DMS decode / encode ===');
  final d = DMS.Decode("48d08'02.7\"N");
  print('Decoded value  = ${d.val}°');                        // 48.134083333333336°
  final enc = DMS.Encode(d.val, DMS.SECOND, 1, DMS.LATITUDE);
  print('Re-encoded DMS = $enc');                             // 48°08'02.7"N
  print('');

  // ── 3. UTM (Karney, ≤ 5 nm accuracy) ──────────────────────────────────────
  print('=== UTM Karney (Berlin 52.52°N 13.41°E) ===');
  final utm = KarneyUTMConverter.fromLatLon(52.52, 13.41);
  print('Zone:      ${utm.zone}${utm.hemisphere}');           // 33N
  print('Easting:   ${utm.easting.toStringAsFixed(3)} m');   // 392118.487 m
  print('Northing:  ${utm.northing.toStringAsFixed(3)} m');  // 5820064.675 m
  print('γ (conv):  ${utm.gamma.toStringAsFixed(9)}°');      // -1.261890431°
  print('k (scale): ${utm.k.toStringAsFixed(12)}');          // 0.999742858842
  print('');

  // ── 4. Transverse Mercator forward (Krüger series) ─────────────────────────
  // Central meridian 9° (= UTM zone 32), point: lat=52°, lon=13.41°
  print('=== KruegerTM forward (lon0=9°, lat=52°, lon=13.41°) ===');
  final tm = KruegerTM(); // WGS84 defaults, k0 = 1.0
  final fwd = tm.forward(9.0, 52.0, 13.41);
  print('x (easting from CM):  ${fwd.x.toStringAsFixed(6)} m');    // 302797.714552 m
  print('y (northing):         ${fwd.y.toStringAsFixed(6)} m');     // 5772534.263476 m
  print('γ (convergence):      ${fwd.gamma.toStringAsFixed(9)}°');  // 3.477748932°
  print('k (scale):            ${fwd.k.toStringAsFixed(12)}');      // 1.001125279691
}
