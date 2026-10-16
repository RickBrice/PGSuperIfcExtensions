# IFC import validation

- Run: 2026-10-16 07:55
- Configuration: Regression:Regression

| Run | Result | Match | Default match | Mismatch | Not imported | Not comparable | Missing |
|---|---|---|---|---|---|---|---|
| [PGSuper-Display](PGSuper-Display.md) | 4 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-Display.usBridge-ids.txt) | 130 of 253 | 115 | 0 | 8 | 0 | 0 |
| [PGSuper-System](PGSuper-System.md) | 4 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-System.usBridge-ids.txt) | 130 of 253 | 115 | 0 | 8 | 0 | 0 |
| [PGSuper-Skew-System](PGSuper-Skew-System.md) | 5 s, export same as baseline, IDS same as baseline, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-System.usBridge-ids.txt) | 30 of 312 | 238 | 8 | 36 | 0 | 0 |
| [PGSuper-Skew-Bearings-System](PGSuper-Skew-Bearings-System.md) | 9 s, export same as baseline, IDS same as baseline, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-Bearings-System.usBridge-ids.txt) | 58 of 312 | 230 | 0 | 24 | 0 | 0 |
| [PGSuper-Skew-AlongPier-System](PGSuper-Skew-AlongPier-System.md) | 6 s, export same as baseline, IDS same as baseline, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-AlongPier-System.usBridge-ids.txt) | 59 of 312 | 229 | 0 | 24 | 0 | 0 |
| [PGSuper-Superelevation-System](PGSuper-Superelevation-System.md) | 4 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-Superelevation-System.usBridge-ids.txt) | 143 of 253 | 102 | 0 | 8 | 0 | 0 |
| [PGSuper-Superelevation-NoProfile-System](PGSuper-Superelevation-NoProfile-System.md) | 4 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-Superelevation-NoProfile-System.usBridge-ids.txt), removed IfcAlignmentVertical before the import | 143 of 253 | 102 | 0 | 8 | 0 | 0 |
| [PennDOT](PennDOT.md) | 10 s, no expected values: 68 of 108 values set by the import | | | | | | |
| [PennDOT-StandardTable](PennDOT-StandardTable.md) | 9 s, log check: 2 of 2 found, no expected values: 58 of 108 values set by the import | | | | | | |
| [Iowa](Iowa.md) | 14 s, no expected values: 208 of 274 values set by the import | | | | | | |
| [Iowa-StandardTable](Iowa-StandardTable.md) | 13 s, log check: 3 of 3 found, no expected values: 146 of 274 values set by the import | | | | | | |
| Standard table <-> IDS | OK: same: 326 declarations (material property sets not compared) | | | | | | |
| Standard table vs usBridge workbook | DIFFERS: 4 CLASS, 3 REQUIRED, 38 NO-VALUE, 13 VALUE, 54 TYPE, 1 NAME, 17 ENUM (20 PSET, 23 RESTRICTED, 36 OPTIONAL, 27 EXTRA, 1 NO-CLASS) (see Standard-vs-usBridge-workbook.txt) | | | | | | |
