# Commit And Push Checklist

## Before Commit

Run from:

```powershell
cd D:\works\work-flutter\geographiclib
```

Required checks:

```powershell
dart test test\triaxial_test.dart
dart run tool\geod3test.dart --quick
```

Optional stronger checks:

```powershell
dart run tool\geod3test.dart --start 5728 --count 1 --inv
dart run tool\geod3test.dart --start 7928 --count 3 --dir
dart run tool\geod3test.dart --start 0 --count 50000 --report-every 5000
```

## Review Scope

Expected transferred areas:

- `lib/triaxial.dart`
- `lib/src/triaxial/`
- `test/triaxial_test.dart`
- `test/data/Geod3Test-v1.txt`
- `tool/geod3test.dart`
- `README_triaxial.md`
- `README_triaxial_transfer.md`
- `CHANGELOG.md`
- `GITHUB_PR_DESCRIPTION.md`

## Suggested Commit Message

```text
Add triaxial ellipsoid support for 2.0.0
```

## Suggested Push Flow

```powershell
git status --short
git add lib/triaxial.dart lib/src/triaxial test/triaxial_test.dart test/data/Geod3Test-v1.txt tool/geod3test.dart README_triaxial.md README_triaxial_transfer.md CHANGELOG.md GITHUB_PR_DESCRIPTION.md COMMIT_PUSH_CHECKLIST.md
git commit -m "Add triaxial ellipsoid support for 2.0.0"
git push
```
