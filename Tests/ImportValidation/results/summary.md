# IFC import validation

- Run: 2026-10-15 14:35
- Configuration: Regression:Regression

| Run | Result | Match | Default match | Mismatch | Not imported | Not comparable | Missing |
|---|---|---|---|---|---|---|---|
| [PGSuper-Display](PGSuper-Display.md) | 4 s, export baseline saved, IDS baseline saved, usBridge IDS: 9 of 599 specifications fail (see PGSuper-Display.usBridge-ids.txt) | 124 of 239 | 101 | 0 | 14 | 0 | 0 |
| [PGSuper-System](PGSuper-System.md) | 8 s, export baseline saved, IDS baseline saved, usBridge IDS: 9 of 599 specifications fail (see PGSuper-System.usBridge-ids.txt) | 124 of 239 | 101 | 0 | 14 | 0 | 0 |
| [PGSuper-Skew-System](PGSuper-Skew-System.md) | 5 s, export baseline saved, IDS baseline saved, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-System.usBridge-ids.txt) | 30 of 294 | 220 | 8 | 36 | 0 | 0 |
| [PGSuper-Skew-Bearings-System](PGSuper-Skew-Bearings-System.md) | 5 s, export baseline saved, IDS baseline saved, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-Bearings-System.usBridge-ids.txt) | 58 of 294 | 212 | 0 | 24 | 0 | 0 |
| [PGSuper-Skew-AlongPier-System](PGSuper-Skew-AlongPier-System.md) | 6 s, export baseline saved, IDS baseline saved, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-AlongPier-System.usBridge-ids.txt) | 59 of 294 | 211 | 0 | 24 | 0 | 0 |
| [PennDOT](PennDOT.md) | 2 s, no expected values: 59 of 102 values set by the import | | | | | | |
| [PennDOT-StandardTable](PennDOT-StandardTable.md) | 2 s, log check: 2 of 2 found, no expected values: 45 of 102 values set by the import | | | | | | |
| [Iowa](Iowa.md) | 8 s, no expected values: 178 of 260 values set by the import | | | | | | |
| [Iowa-StandardTable](Iowa-StandardTable.md) | 7 s, log check: 3 of 3 found, no expected values: 118 of 260 values set by the import | | | | | | |
| Standard table <-> IDS | OK: same: 326 declarations (material property sets not compared) | | | | | | |
| Standard table vs usBridge workbook | DIFFERS: 4 CLASS, 3 REQUIRED, 38 NO-VALUE, 13 VALUE, 54 TYPE, 1 NAME, 17 ENUM (20 PSET, 23 RESTRICTED, 36 OPTIONAL, 27 EXTRA, 1 NO-CLASS) (see Standard-vs-usBridge-workbook.txt) | | | | | | |
