/// UTM ↔ geographic conversions.
///
/// Two implementations:
/// - [KarneyUTM] — Karney/Krüger 6th-order series, ≤ 5 nm accuracy.
/// - [UTMConverter] — Snyder formulas, ~1 mm accuracy.
///
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

library utm;

export 'src/krueger_tm.dart';
export 'src/karney_utm.dart';
export 'src/utm.dart';
