# IFC import validation

- Run: 2026-10-15 12:50
- Configuration: Regression:Regression

| Run | Result | Match | Default match | Mismatch | Not imported | Not comparable | Missing |
|---|---|---|---|---|---|---|---|
| [PGSuper-Display](PGSuper-Display.md) | 3 s, export same as baseline, IDS same as baseline, usBridge IDS: 15 of 599 specifications fail (see PGSuper-Display.usBridge-ids.txt) | 124 of 239 | 101 | 0 | 14 | 0 | 0 |
| [PGSuper-System](PGSuper-System.md) | 3 s, export same as baseline, IDS same as baseline, usBridge IDS: 15 of 599 specifications fail (see PGSuper-System.usBridge-ids.txt) | 124 of 239 | 101 | 0 | 14 | 0 | 0 |
| [PGSuper-Skew-System](PGSuper-Skew-System.md) | 5 s, export same as baseline, IDS same as baseline, usBridge IDS: 11 of 599 specifications fail (see PGSuper-Skew-System.usBridge-ids.txt) | 30 of 294 | 220 | 8 | 36 | 0 | 0 |
| [PGSuper-Skew-Bearings-System](PGSuper-Skew-Bearings-System.md) | 5 s, export same as baseline, IDS same as baseline, usBridge IDS: 12 of 599 specifications fail (see PGSuper-Skew-Bearings-System.usBridge-ids.txt) | 58 of 294 | 212 | 0 | 24 | 0 | 0 |
| [PGSuper-Skew-AlongPier-System](PGSuper-Skew-AlongPier-System.md) | 5 s, export same as baseline, IDS same as baseline, usBridge IDS: 12 of 599 specifications fail (see PGSuper-Skew-AlongPier-System.usBridge-ids.txt) | 59 of 294 | 211 | 0 | 24 | 0 | 0 |
| [PennDOT](PennDOT.md) | 2 s, no expected values: 59 of 102 values set by the import | | | | | | |
| [PennDOT-StandardTable](PennDOT-StandardTable.md) | 2 s, log check: 2 of 2 found, no expected values: 45 of 102 values set by the import | | | | | | |
| [Iowa](Iowa.md) | 7 s, no expected values: 178 of 260 values set by the import | | | | | | |
| [Iowa-StandardTable](Iowa-StandardTable.md) | 7 s, log check: 3 of 3 found, no expected values: 118 of 260 values set by the import | | | | | | |
| Standard table <-> IDS | OK: same: 273 declarations (material property sets not compared) | | | | | | |
| Standard table vs usBridge workbook | DIFFERS: 4 CLASS, 3 REQUIRED, 38 NO-VALUE, 13 VALUE, 52 TYPE, 1 NAME, 17 ENUM (23 PSET, 23 RESTRICTED, 22 OPTIONAL, 26 EXTRA, 2 NO-CLASS) (see Standard-vs-usBridge-workbook.txt) | | | | | | |
