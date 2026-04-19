// dms.dart
// Dart port of DMS.js from geographiclib-dms v2.2.0
//
// ==========================================================================
// Original C++/JS: Copyright (c) Charles Karney (2008-2024)
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

// Hemisphere / format indicator values (mirrors DMS.NONE etc. in JS).
// Kept as a plain class with static constants so callers can use
// DMS.NONE, DMS.LATITUDE, etc. — identical naming to the JS API.

/// Dart port of Charles Karney's DMS class from geographiclib-js.
///
/// Provides `Decode`, `Encode`, `DecodeLatLon`, `DecodeAngle`,
/// `DecodeAzimuth` — a faithful port of the JavaScript implementation.
///
/// Also retains the original helper API (`fromDecimal`, `toDecimal`,
/// `format`, `parse`) for backward compatibility.
class DMS {
  DMS._(); // static-only class

  // -------------------------------------------------------------------------
  // Constants
  // -------------------------------------------------------------------------

  /// No hemisphere designator; angle formatted as plain signed number.
  static const int NONE = 0;

  /// N/S hemisphere; format as latitude (2-digit degree pad).
  static const int LATITUDE = 1;

  /// E/W hemisphere; format as longitude (3-digit degree pad).
  static const int LONGITUDE = 2;

  /// Format as azimuth: range [0, 360), 3-digit degree pad, no sign.
  static const int AZIMUTH = 3;

  /// Trailing component is degrees.
  static const int DEGREE = 0;

  /// Trailing component is minutes.
  static const int MINUTE = 1;

  /// Trailing component is seconds.
  static const int SECOND = 2;

  // Internal helpers --------------------------------------------------------

  static const String _hemispheres = 'SNWE';
  static const String _signs = '-+';
  static const String _digits = '0123456789';
  static const String _dmsindicators = "D'\":";
  static const String _dmsindicatorsu = "°'\""; // degree symbol variants

  static int _lookup(String s, String c) =>
      s.indexOf(c.toUpperCase());

  static String _zerofill(String s, int n) {
    final int pad = n - s.length;
    if (pad <= 0) return s;
    return ('0000'.substring(0, pad.clamp(0, 4))) + s;
  }

  // -------------------------------------------------------------------------
  // Decode
  // -------------------------------------------------------------------------

  /// Decodes a DMS string into `{val, ind}`.
  ///
  /// Returns a record `({double val, int ind})` where `ind` is one of
  /// [NONE], [LATITUDE], or [LONGITUDE].
  ///
  /// Accepts all formats supported by the C++/JS library, including
  /// `48d30'40.5"S`, `48:30:40.5`, `127:54:3.123123W`, etc.
  /// Also accepts `nan`, `inf`.
  ///
  /// Throws [FormatException] on illegal input.
  static ({double val, int ind}) Decode(String dms) {
    // Normalise Unicode and alternate symbols to ASCII d / ' / " / + / -
    String dmsa = dms
        .replaceAll('\u00b0', 'd')
        .replaceAll('\u00ba', 'd')
        .replaceAll('\u2070', 'd')
        .replaceAll('\u02da', 'd')
        .replaceAll('\u2218', 'd')
        .replaceAll('*', 'd')
        .replaceAll('`', 'd')
        .replaceAll('\u2032', "'")
        .replaceAll('\u2035', "'")
        .replaceAll('\u00b4', "'")
        .replaceAll('\u2018', "'")
        .replaceAll('\u2019', "'")
        .replaceAll('\u201b', "'")
        .replaceAll('\u02b9', "'")
        .replaceAll('\u02ca', "'")
        .replaceAll('\u02cb', "'")
        .replaceAll('\u2033', '"')
        .replaceAll('\u2036', '"')
        .replaceAll('\u02dd', '"')
        .replaceAll('\u201c', '"')
        .replaceAll('\u201d', '"')
        .replaceAll('\u201f', '"')
        .replaceAll('\u02ba', '"')
        .replaceAll('\u2795', '+')
        .replaceAll('\u2064', '+')
        .replaceAll('\u2010', '-')
        .replaceAll('\u2011', '-')
        .replaceAll('\u2013', '-')
        .replaceAll('\u2014', '-')
        .replaceAll('\u2212', '-')
        .replaceAll('\u2796', '-')
        .replaceAll('\u00a0', '')
        .replaceAll('\u2007', '')
        .replaceAll('\u2009', '')
        .replaceAll('\u200a', '')
        .replaceAll('\u200b', '')
        .replaceAll('\u202f', '')
        .replaceAll('\u2063', '')
        .replaceAll("''", '"')
        .trim();

    final int end = dmsa.length;
    double v = -0.0;
    int i = 0;
    int ind1 = NONE;

    for (int p = 0, pb = 0; p < end; p = pb, ++i) {
      int pa = p;
      // Skip initial hemisphere letter (only for the first piece)
      if (i == 0 && _lookup(_hemispheres, dmsa[pa]) >= 0) ++pa;
      // Skip over initial sign
      if (i > 0 || (pa < end && _lookup(_signs, dmsa[pa]) >= 0)) ++pa;
      // Find next sign
      int mi = dmsa.indexOf('-', pa);
      int pi = dmsa.indexOf('+', pa);
      if (mi < 0) mi = end;
      if (pi < 0) pi = end;
      pb = mi < pi ? mi : pi;

      final decoded = _internalDecode(dmsa.substring(p, pb));
      v += decoded.val;
      final int ind2 = decoded.ind;
      if (ind1 == NONE) {
        ind1 = ind2;
      } else if (ind2 != NONE && ind1 != ind2) {
        throw FormatException(
            'Incompatible hemisphere designators in ${dmsa.substring(0, pb)}');
      }
    }
    if (i == 0) throw FormatException('Empty or incomplete DMS string $dms');
    return (val: v, ind: ind1);
  }

  static ({double val, int ind}) _internalDecode(String dmsa) {
    // Try nan / inf first
    final double? special = _numMatch(dmsa);
    if (special != null) return (val: special, ind: NONE);

    return _parseComponents(dmsa);
  }

  // Separated so we can use early-return instead of do..while(false)
  static ({double val, int ind}) _parseComponents(String dmsa) {
    int ind1 = NONE;
    int sign = 1;
    int beg = 0;
    int end = dmsa.length;
    int k;

    // Leading hemisphere
    if (end > beg && (k = _lookup(_hemispheres, dmsa[beg])) >= 0) {
      ind1 = (k & 2) != 0 ? LONGITUDE : LATITUDE;
      sign = (k & 1) != 0 ? 1 : -1;
      ++beg;
    }
    // Trailing hemisphere
    if (end > beg && (k = _lookup(_hemispheres, dmsa[end - 1])) >= 0) {
      if (ind1 != NONE) {
        if (dmsa[beg - 1].toUpperCase() == dmsa[end - 1].toUpperCase()) {
          throw FormatException(
              'Repeated hemisphere indicators ${dmsa[beg - 1]} in ${dmsa.substring(beg - 1, end)}');
        } else {
          throw FormatException(
              'Contradictory hemisphere indicators ${dmsa[beg - 1]} and ${dmsa[end - 1]} in ${dmsa.substring(beg - 1, end)}');
        }
      }
      ind1 = (k & 2) != 0 ? LONGITUDE : LATITUDE;
      sign = (k & 1) != 0 ? 1 : -1;
      --end;
    }
    // Leading sign
    if (end > beg && (k = _lookup(_signs, dmsa[beg])) >= 0) {
      sign *= (k != 0) ? 1 : -1;
      ++beg;
    }

    if (end == beg) {
      throw FormatException('Empty or incomplete DMS string $dmsa');
    }

    final List<int> ipieces = [0, 0, 0];
    final List<double> fpieces = [0, 0, 0];
    int npiece = 0;
    int icurrent = 0;
    double fcurrent = 0;
    int ncurrent = 0;
    int p = beg;
    bool pointseen = false;
    int digcount = 0;
    int intcount = 0;
    String errormsg = '';

    outer:
    while (p < end) {
      final String x = dmsa[p++];
      int kd;
      if ((kd = _lookup(_digits, x)) >= 0) {
        ++ncurrent;
        if (digcount > 0) {
          ++digcount;
        } else {
          icurrent = 10 * icurrent + kd;
          ++intcount;
        }
      } else if (x == '.') {
        if (pointseen) {
          errormsg = 'Multiple decimal points in ${dmsa.substring(beg, end)}';
          break;
        }
        pointseen = true;
        digcount = 1;
      } else if ((kd = _lookup(_dmsindicators, x)) >= 0) {
        if (kd >= 3) {
          if (p == end) {
            errormsg =
                'Illegal for colon to appear at the end of ${dmsa.substring(beg, end)}';
            break outer;
          }
          kd = npiece;
        }
        if (kd == npiece - 1) {
          errormsg =
              'Repeated ${_componentName(kd)} component in ${dmsa.substring(beg, end)}';
          break;
        } else if (kd < npiece) {
          errormsg =
              '${_componentName(kd)} component follows ${_componentName(npiece - 1)} component in ${dmsa.substring(beg, end)}';
          break;
        }
        if (ncurrent == 0) {
          errormsg =
              'Missing numbers in ${_componentName(kd)} component of ${dmsa.substring(beg, end)}';
          break;
        }
        if (digcount > 0) {
          fcurrent = double.parse(
              dmsa.substring(p - intcount - digcount - 1, p - 1));
          icurrent = 0;
        }
        ipieces[kd] = icurrent;
        fpieces[kd] = icurrent + fcurrent;
        if (p < end) {
          npiece = kd + 1;
          if (npiece >= 3) {
            errormsg =
                'More than 3 DMS components in ${dmsa.substring(beg, end)}';
            break;
          }
          icurrent = 0;
          fcurrent = 0;
          ncurrent = 0;
          digcount = 0;
          intcount = 0;
        }
      } else if (_lookup(_signs, x) >= 0) {
        errormsg = 'Internal sign in DMS string ${dmsa.substring(beg, end)}';
        break;
      } else {
        errormsg =
            'Illegal character $x in DMS string ${dmsa.substring(beg, end)}';
        break;
      }
    }

    if (errormsg.isEmpty && _lookup(_dmsindicators, dmsa[p - 1]) < 0) {
      if (npiece >= 3) {
        errormsg =
            'Extra text following seconds in DMS string ${dmsa.substring(beg, end)}';
      } else if (ncurrent == 0) {
        errormsg =
            'Missing numbers in trailing component of ${dmsa.substring(beg, end)}';
      } else {
        if (digcount > 0) {
          fcurrent =
              double.parse(dmsa.substring(p - intcount - digcount, p));
          icurrent = 0;
        }
        ipieces[npiece] = icurrent;
        fpieces[npiece] = icurrent + fcurrent;
      }
    }

    if (errormsg.isEmpty && pointseen && digcount == 0) {
      errormsg =
          'Decimal point in non-terminal component of ${dmsa.substring(beg, end)}';
    }
    if (errormsg.isEmpty && (ipieces[1] >= 60 || fpieces[1] > 60)) {
      errormsg = 'Minutes ${fpieces[1]} not in range [0,60)';
    }
    if (errormsg.isEmpty && (ipieces[2] >= 60 || fpieces[2] > 60)) {
      errormsg = 'Seconds ${fpieces[2]} not in range [0,60)';
    }

    if (errormsg.isNotEmpty) throw FormatException(errormsg);

    final double computed = sign *
        (fpieces[2] != 0
            ? (60 * (60 * fpieces[0] + fpieces[1]) + fpieces[2]) / 3600
            : (fpieces[1] != 0
                ? (60 * fpieces[0] + fpieces[1]) / 60
                : fpieces[0]));
    return (val: computed, ind: ind1);
  }

  static String _componentName(int k) =>
      const ['degrees', 'minutes', 'seconds'][k];

  static double? _numMatch(String s) {
    if (s.length < 3) return null;
    String t = s.toUpperCase().replaceAll(RegExp(r'0+$'), '');
    final int sign = t.startsWith('-') ? -1 : 1;
    final int p0 = (t.startsWith('-') || t.startsWith('+')) ? 1 : 0;
    if (t.length - p0 < 3) return null;
    t = t.substring(p0);
    if (t == 'NAN' ||
        t == '1.#QNAN' ||
        t == '1.#SNAN' ||
        t == '1.#IND' ||
        t == '1.#R') {
      return double.nan;
    } else if (t == 'INF' || t == '1.#INF' || t == 'INFINITY') {
      return sign * double.infinity;
    }
    return null;
  }

  // -------------------------------------------------------------------------
  // DecodeLatLon
  // -------------------------------------------------------------------------

  /// Decodes two DMS strings as a latitude/longitude pair.
  ///
  /// Returns `({double lat, double lon})`.
  /// If [longfirst] is true, the first string is treated as longitude
  /// when neither string has a hemisphere designator.
  static ({double lat, double lon}) DecodeLatLon(String stra, String strb,
      {bool longfirst = false}) {
    final va = Decode(stra);
    final vb = Decode(strb);
    double a = va.val;
    int ia = va.ind;
    double b = vb.val;
    int ib = vb.ind;

    if (ia == NONE && ib == NONE) {
      ia = longfirst ? LONGITUDE : LATITUDE;
      ib = longfirst ? LATITUDE : LONGITUDE;
    } else if (ia == NONE) {
      ia = LATITUDE + LONGITUDE - ib;
    } else if (ib == NONE) {
      ib = LATITUDE + LONGITUDE - ia;
    }
    if (ia == ib) {
      throw FormatException(
          'Both $stra and $strb interpreted as ${ia == LATITUDE ? "latitudes" : "longitudes"}');
    }
    final double lat = ia == LATITUDE ? a : b;
    final double lon = ia == LATITUDE ? b : a;
    if (lat.abs() > 90) {
      throw FormatException('Latitude $lat not in [-90,90]');
    }
    return (lat: lat, lon: lon);
  }

  // -------------------------------------------------------------------------
  // DecodeAngle
  // -------------------------------------------------------------------------

  /// Decodes a DMS string as an arc angle (no hemisphere allowed).
  static double DecodeAngle(String angstr) {
    final r = Decode(angstr);
    if (r.ind != NONE) {
      throw FormatException(
          'Arc angle $angstr includes a hemisphere N/E/W/S');
    }
    return r.val;
  }

  // -------------------------------------------------------------------------
  // DecodeAzimuth
  // -------------------------------------------------------------------------

  /// Decodes a DMS string as an azimuth (E/W allowed, N/S not allowed).
  static double DecodeAzimuth(String azistr) {
    final r = Decode(azistr);
    if (r.ind == LATITUDE) {
      throw FormatException('Azimuth $azistr has a latitude hemisphere N/S');
    }
    return r.val;
  }

  // -------------------------------------------------------------------------
  // Encode
  // -------------------------------------------------------------------------

  /// Encodes an angle (degrees) to a DMS string.
  ///
  /// [trailing] is one of [DEGREE], [MINUTE], or [SECOND].
  /// [prec] is the number of decimal digits for the trailing component.
  /// [ind] is one of [NONE], [LATITUDE], [LONGITUDE], or [AZIMUTH].
  /// [dmssep] if non-null, replaces the `°`, `'`, `"` separators.
  static String Encode(double angle, int trailing, int prec,
      [int ind = NONE, String? dmssep]) {
    if (!angle.isFinite) {
      if (angle.isNaN) return 'nan';
      return angle < 0 ? '-inf' : 'inf';
    }

    if (angle.abs() >= 1e21) {
      // toFixed not reliable beyond 1e21 — format in exponential
      return angle.toStringAsExponential().replaceAll('e+', 'e');
    }

    final bool usesep = dmssep != null;
    final String sep = dmssep ?? '';

    if (ind == AZIMUTH) {
      angle %= 360;
      if (angle < 0) angle += 360;
      // convert -0 to +0
      if (angle == 0) angle = 0.0;
    }

    final int sign =
        (angle < 0 || (angle == 0 && angle.isNegative)) ? -1 : 1;
    angle *= sign;

    // precision cap: 15 - 2*trailing gives full real precision in [-90,90]
    final int effectivePrec = prec.clamp(0, 15 - 2 * trailing);

    double scale = 1.0;
    for (int i = 0; i < trailing; ++i) scale *= 60;

    final int idegree = trailing == DEGREE ? 0 : angle.floor();
    double fdegree = (angle - idegree) * scale;

    // Format fdegree with effectivePrec decimal places
    String s = fdegree.toStringAsFixed(effectivePrec);

    String degree, minute = '', second = '';

    switch (trailing) {
      case DEGREE:
        degree = s;
        break;
      default:
        // MINUTE or SECOND
        int dotPos = s.indexOf('.');
        int ipart;
        String fpart;
        if (dotPos < 0) {
          ipart = int.parse(s);
          fpart = '';
        } else if (dotPos == 0) {
          ipart = 0;
          fpart = s;
        } else {
          ipart = int.parse(s.substring(0, dotPos));
          fpart = s.substring(dotPos);
        }
        if (trailing == MINUTE) {
          minute = '${ipart % 60}$fpart';
          degree = (ipart ~/ 60 + idegree).toStringAsFixed(0);
        } else {
          // SECOND
          second = '${ipart % 60}$fpart';
          int tmp = ipart ~/ 60;
          minute = '${tmp % 60}';
          degree = (tmp ~/ 60 + idegree).toStringAsFixed(0);
        }
        break;
    }

    // Extra width for decimal point
    final int dp = effectivePrec > 0 ? effectivePrec + 1 : 0;

    // degree pad width: NONE→0, LATITUDE→2, LONGITUDE→3, AZIMUTH→3
    int degWidth = ind == NONE ? 0 : (ind == LATITUDE ? 2 : 3);

    StringBuffer buf = StringBuffer();
    if (ind == NONE && sign < 0) buf.write('-');

    switch (trailing) {
      case DEGREE:
        buf.write(_zerofill(degree, degWidth + dp));
        if (!usesep) buf.write(_dmsindicatorsu[0]);
        break;
      case MINUTE:
        buf.write(_zerofill(degree, degWidth));
        buf.write(usesep ? sep : _dmsindicatorsu[0]);
        buf.write(_zerofill(minute, 2 + dp));
        if (!usesep) buf.write(_dmsindicatorsu[1]);
        break;
      default: // SECOND
        buf.write(_zerofill(degree, degWidth));
        buf.write(usesep ? sep : _dmsindicatorsu[0]);
        buf.write(_zerofill(minute, 2));
        buf.write(usesep ? sep : _dmsindicatorsu[1]);
        buf.write(_zerofill(second, 2 + dp));
        if (!usesep) buf.write(_dmsindicatorsu[2]);
        break;
    }

    if (ind != NONE && ind != AZIMUTH) {
      // hemispheres = 'SNWE'
      // LATITUDE: sign<0 → 'S'(0), sign>0 → 'N'(1)
      // LONGITUDE: sign<0 → 'W'(2), sign>0 → 'E'(3)
      final int hi = (ind == LATITUDE ? 0 : 2) + (sign < 0 ? 0 : 1);
      buf.write(_hemispheres[hi]);
    }

    return buf.toString();
  }

  // =========================================================================
  // Legacy / convenience API (backward-compatible)
  // =========================================================================

  // -------------------------------------------------------------------------
  // Decimal → DMS record
  // -------------------------------------------------------------------------

  /// Converts [decimal] degrees to a (degrees, minutes, seconds) record.
  static ({int degrees, int minutes, double seconds}) fromDecimal(
      double decimal) {
    final bool negative = decimal < 0;
    final double abs = decimal.abs();
    final int deg = abs.floor().toInt();
    final double minFull = (abs - deg) * 60.0;
    final int min = minFull.floor().toInt();
    final double sec = (minFull - min) * 60.0;
    return (
      degrees: negative ? -deg : deg,
      minutes: min,
      seconds: sec,
    );
  }

  // -------------------------------------------------------------------------
  // DMS → Decimal
  // -------------------------------------------------------------------------

  /// Converts degrees, minutes, seconds and an optional hemisphere to decimal.
  static double toDecimal(
    int degrees,
    int minutes,
    double seconds, [
    String hemisphere = 'N',
  ]) {
    final double decimal = degrees.abs() + minutes / 60.0 + seconds / 3600.0;
    final String h = hemisphere.toUpperCase();
    if (h == 'S' || h == 'W' || degrees < 0) return -decimal;
    return decimal;
  }

  // -------------------------------------------------------------------------
  // Simple formatter
  // -------------------------------------------------------------------------

  /// Formats [decimal] degrees as a `48°51'23.76"N` style string.
  static String format(double decimal, {bool isLatitude = true}) {
    final bool negative = decimal < 0;
    final dms = fromDecimal(decimal.abs());
    final String dir =
        isLatitude ? (negative ? 'S' : 'N') : (negative ? 'W' : 'E');
    return "${dms.degrees}°${dms.minutes}'${dms.seconds.toStringAsFixed(2)}\"$dir";
  }

  // -------------------------------------------------------------------------
  // Simple parser
  // -------------------------------------------------------------------------

  /// Parses a DMS string to decimal degrees (simple regex-based).
  ///
  /// For the full Karney parser, use [Decode].
  static double parse(String dms) {
    final String s = dms.trim();
    final bool leadingMinus = s.startsWith('-');
    final String abs = leadingMinus ? s.substring(1) : s;

    final RegExp re = RegExp(
        r"(\d+)[°:\s]+(\d+)[':\s]+([\d.]+)[^NSEWnsew0-9]*([NSEWnsew]?)");
    final Match? m = re.firstMatch(abs);
    if (m == null) throw FormatException('Cannot parse DMS string: "$dms"');

    final int deg = int.parse(m.group(1)!);
    final int min = int.parse(m.group(2)!);
    final double sec = double.parse(m.group(3)!);
    final String hem = m.group(4)!.toUpperCase();

    double decimal = deg + min / 60.0 + sec / 3600.0;
    if (leadingMinus || hem == 'S' || hem == 'W') decimal = -decimal;
    return decimal;
  }
}
