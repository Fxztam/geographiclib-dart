import 'dart:io';
import 'dart:math' as math;

import 'package:geographiclib_dart/triaxial.dart';

const double _eps = 2.220446049250313e-16;
const double _invThreshold = 500.0;
const double _invMeanThreshold = 20.0;
const double _dirThreshold = 5000.0;

class _FailureRecord {
  final int lineIdx;
  final String type;
  final Map<String, Object?> details;

  _FailureRecord(this.lineIdx, this.type, this.details);

  @override
  String toString() {
    final parts = <String>['line=$lineIdx', 'type=$type'];
    details.forEach((key, value) => parts.add('$key=$value'));
    return parts.join(' ');
  }
}

class _Stats {
  int n = 0;
  double errsMax = 0;
  double errsSum = 0;
  double errAlp1Max = 0;
  double errAlp2Max = 0;
  double errDirMax = 0;
  double errDirSum = 0;
  final List<_FailureRecord> failures = <_FailureRecord>[];
}

String _fmt(double value, [int digits = 2]) {
  if (!value.isFinite) return value.toString();
  return value.toStringAsFixed(digits);
}

double _angleDiffDeg(double a, double b) {
  final d = (a - b).abs() % 360.0;
  return d > 180.0 ? 360.0 - d : d;
}

double _maxComponentError(List<double> actual, List<double> expected) {
  var maxErr = 0.0;
  for (var index = 0; index < actual.length; index++) {
    maxErr = math.max(maxErr, (actual[index] - expected[index]).abs() / _eps);
  }
  return maxErr;
}

Angle _copyAngle(Angle angle) => angle.copy();

void _recordFailure(
    _Stats stats, int lineIdx, String type, Map<String, Object?> details) {
  if (stats.failures.length < 5) {
    stats.failures.add(_FailureRecord(lineIdx, type, details));
  }
}

void main(List<String> args) {
  var startIndex = 0;
  int? countLimit;
  var doInverse = true;
  var doDirect = true;
  var reportEvery = 10000;
  String? dataPathArg;

  for (var index = 0; index < args.length; index++) {
    final arg = args[index];
    if (arg == '--quick') {
      countLimit = 1000;
    } else if (arg == '--data' && index + 1 < args.length) {
      dataPathArg = args[++index];
    } else if (arg == '--start' && index + 1 < args.length) {
      startIndex = math.max(0, int.parse(args[++index]));
    } else if (arg == '--count' && index + 1 < args.length) {
      countLimit = int.parse(args[++index]);
    } else if (arg == '--report-every' && index + 1 < args.length) {
      reportEvery = math.max(0, int.parse(args[++index]));
    } else if (arg == '--inv') {
      doDirect = false;
    } else if (arg == '--dir') {
      doInverse = false;
    } else {
      stderr.writeln('Unknown argument: $arg');
      exitCode = 2;
      return;
    }
  }

  final envDataPath = Platform.environment['GEOD3TEST_DATA'];
  final dataCandidates = <File>[
    if (dataPathArg != null) File(dataPathArg),
    if (envDataPath != null && envDataPath.isNotEmpty) File(envDataPath),
    File.fromUri(Platform.script.resolve('../test/data/Geod3Test-v1.txt')),
    File.fromUri(Platform.script.resolve('../../test/data/Geod3Test-v1.txt')),
  ];
  final dataFile = dataCandidates.firstWhere(
    (file) => file.existsSync(),
    orElse: () => dataCandidates.first,
  );
  if (!dataFile.existsSync()) {
    stderr.writeln('Reference dataset not found.');
    stderr.writeln('Tried: ${dataFile.path}');
    stderr.writeln(
        'Provide it via --data <path> or GEOD3TEST_DATA, or place Geod3Test-v1.txt under test\\data\\.');
    stderr.writeln(
        'Origin: https://doi.org/10.5281/zenodo.12510796  (Zenodo dataset, version 1.0)');
    exitCode = 2;
    return;
  }

  final lines = dataFile.readAsLinesSync();
  final total =
      math.min(lines.length, startIndex + (countLimit ?? lines.length));
  final ell = Ellipsoid3(math.sqrt(2.0), 1.0, 1.0 / math.sqrt(2.0));
  final geod = Geodesic3(ell);
  final stats = _Stats();

  stdout.writeln(
      'Ellipsoid: a=${_fmt(math.sqrt(2.0), 8)} b=1 c=${_fmt(1.0 / math.sqrt(2.0), 8)}');
  stdout.writeln('k2=${_fmt(geod.k2, 8)} kp2=${_fmt(geod.kp2, 8)}');

  final startedAt = DateTime.now();
  for (var lineNumber = startIndex; lineNumber < total; lineNumber++) {
    final line = lines[lineNumber].trim();
    if (line.isEmpty || line.startsWith('#')) continue;

    final cols = line.split(RegExp(r'\s+'));
    if (cols.length < 7) continue;

    final bet1r = double.parse(cols[0]);
    final omg1r = double.parse(cols[1]);
    final alp1Ref = double.parse(cols[2]);
    final bet2r = double.parse(cols[3]);
    final omg2r = double.parse(cols[4]);
    final alp2Ref = double.parse(cols[5]);
    final s12Ref = double.parse(cols[6]);

    stats.n++;

    if (doInverse) {
      try {
        final result = geod.inverse(
          Angle.fromDegrees(bet1r),
          Angle.fromDegrees(omg1r),
          Angle.fromDegrees(bet2r),
          Angle.fromDegrees(omg2r),
        );
        final errs = (result.s12 - s12Ref).abs() / _eps;
        stats.errsMax = math.max(stats.errsMax, errs);
        stats.errsSum += errs;

        final ea1 = _angleDiffDeg(result.alp1.degrees0(), alp1Ref);
        final ea2 = _angleDiffDeg(result.alp2.degrees0(), alp2Ref);
        stats.errAlp1Max = math.max(stats.errAlp1Max, ea1);
        stats.errAlp2Max = math.max(stats.errAlp2Max, ea2);

        if (errs > _invThreshold) {
          _recordFailure(stats, lineNumber + 1, 'inv_dist', {
            'bet1': bet1r,
            'omg1': omg1r,
            'bet2': bet2r,
            'omg2': omg2r,
            's12': result.s12,
            's12_ref': s12Ref,
            'errs_eps': errs,
            'alp1': result.alp1.degrees0(),
            'alp1_ref': alp1Ref,
            'alp2': result.alp2.degrees0(),
            'alp2_ref': alp2Ref,
          });
        }
      } catch (error) {
        _recordFailure(stats, lineNumber + 1, 'inv_exception', {
          'message': error.toString(),
        });
      }
    }

    if (doDirect) {
      try {
        final line3 = GeodesicLine3.fromAngles(
          geod,
          Angle.fromDegrees(bet1r),
          Angle.fromDegrees(omg1r),
          Angle.fromDegrees(alp1Ref),
        );
        final pos = line3.position(s12Ref);
        final refBet = Angle.fromDegrees(bet2r);
        final refOmg = Angle.fromDegrees(omg2r);
        final candidates = <({Angle bet2, Angle omg2, Angle alp2})>[];

        final p0 = (
          bet2: _copyAngle(pos.bet2),
          omg2: _copyAngle(pos.omg2),
          alp2: _copyAngle(pos.alp2)
        );
        Ellipsoid3.angNorm(p0.bet2, p0.omg2, p0.alp2);
        candidates.add(p0);

        if ((bet2r.abs() - 90.0).abs() < 1e-12 ||
            (p0.bet2.degrees0().abs() - 90.0).abs() < 1e-12) {
          final p1 = (
            bet2: _copyAngle(pos.bet2),
            omg2: _copyAngle(pos.omg2),
            alp2: _copyAngle(pos.alp2)
          );
          Ellipsoid3.angNorm(p1.bet2, p1.omg2, p1.alp2, true);
          candidates.add(p1);
        }

        var errDir = double.infinity;
        for (final candidate in candidates) {
          final got = ell.elliptocart2(candidate.bet2, candidate.omg2);
          final ref = ell.elliptocart2(refBet, refOmg);
          errDir = math.min(errDir, _maxComponentError(got, ref));
        }
        stats.errDirMax = math.max(stats.errDirMax, errDir);
        stats.errDirSum += errDir;

        if (errDir > _dirThreshold) {
          _recordFailure(stats, lineNumber + 1, 'dir_pos', {
            'bet1': bet1r,
            'omg1': omg1r,
            'bet2': bet2r,
            'omg2': omg2r,
            'alp1_ref': alp1Ref,
            'alp2_ref': alp2Ref,
            'errdir_eps': errDir,
            'pos_bet2': pos.bet2.degrees0(),
            'pos_omg2': pos.omg2.degrees0(),
            'pos_alp2': pos.alp2.degrees0(),
          });
        }
      } catch (error) {
        _recordFailure(stats, lineNumber + 1, 'dir_exception', {
          'message': error.toString(),
        });
      }
    }

    if (reportEvery > 0 && ((lineNumber - startIndex + 1) % reportEvery) == 0) {
      stdout.write(
        '  ${lineNumber - startIndex + 1}/${total - startIndex}'
        '  inv_s12_max=${_fmt(stats.errsMax, 1)} eps'
        '  dir_max=${_fmt(stats.errDirMax, 1)} eps\r',
      );
    }
  }

  final elapsed = DateTime.now().difference(startedAt).inMilliseconds / 1000.0;
  stdout.writeln(
      '\n\nResults after ${stats.n} test cases from lines ${startIndex + 1}-$total (${elapsed.toStringAsFixed(1)}s):');
  var ok = true;

  if (doInverse) {
    final invMean = stats.n == 0 ? 0.0 : stats.errsSum / stats.n;
    final invMaxPass = stats.errsMax <= _invThreshold;
    final invMeanPass = invMean <= _invMeanThreshold;
    stdout.writeln('  Inverse distance (errs):');
    stdout.writeln(
        '    max = ${_fmt(stats.errsMax)} eps  (pass: <= ${_fmt(_invThreshold, 0)})  ${invMaxPass ? 'PASS' : 'FAIL'}');
    stdout.writeln(
        '    mean = ${invMean.toStringAsFixed(3)} eps  (pass: <= ${_fmt(_invMeanThreshold, 0)})  ${invMeanPass ? 'PASS' : 'FAIL'}');
    stdout.writeln(
        '  Azimuth errors:  alp1_max=${(stats.errAlp1Max * 1e12).toStringAsFixed(3)} peps   alp2_max=${(stats.errAlp2Max * 1e12).toStringAsFixed(3)} peps');
    ok = ok && invMaxPass && invMeanPass;
  }

  if (doDirect) {
    final dirPass = stats.errDirMax <= _dirThreshold;
    stdout.writeln('  Direct position (errdir):');
    stdout.writeln(
        '    max = ${_fmt(stats.errDirMax)} eps  (pass: <= ${_fmt(_dirThreshold, 0)})  ${dirPass ? 'PASS' : 'FAIL'}');
    ok = ok && dirPass;
  }

  if (stats.failures.isNotEmpty) {
    stdout.writeln('\n  First failures:');
    for (final failure in stats.failures) {
      stdout.writeln('    $failure');
    }
  }

  exitCode = ok ? 0 : 1;
}
