// GeographicLib Triaxial Geodesic — Dart port.
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
//
// Exports the full triaxial geodesic API:
// - [Angle] — angle with turn tracking
// - [EllipticFunction3] — Carlson symmetric integrals + Jacobi functions
// - [Trigfun] / [TrigfunExt] — Fourier approximation of integrands
// - [Ellipsoid3] — triaxial ellipsoid parameterization
// - [GeodesicLine3] — triaxial geodesic line (direct problem)
// - [Geodesic3] — triaxial geodesic solver (inverse + direct)
library triaxial;

export 'src/triaxial/angle.dart';
export 'src/triaxial/elliptic_function3.dart';
export 'src/triaxial/trigfun.dart';
export 'src/triaxial/ellipsoid3.dart';
export 'src/triaxial/geodesic_line3.dart';
export 'src/triaxial/geodesic3.dart';
