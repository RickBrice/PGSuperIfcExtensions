# IFC import validation

- Run: 2026-10-20 15:34
- Configuration: Regression

| Run | Result | Match | Default match | Mismatch | Not imported | Not comparable | Missing |
|---|---|---|---|---|---|---|---|
| [PGSuper-Display](PGSuper-Display.md) | 25 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-Display.usBridge-ids.txt) | 306 of 473 | 154 | 13 | 0 | 0 | 0 |
| [PGSuper-System](PGSuper-System.md) | 21 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-System.usBridge-ids.txt) | 306 of 473 | 154 | 13 | 0 | 0 | 0 |
| [PGSuper-WSDOT-System](PGSuper-WSDOT-System.md) | 14 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-WSDOT-System.usBridge-ids.txt), configuration "WSDOT (Local)":"WSDOT (Local)" | 319 of 473 | 154 | 0 | 0 | 0 | 0 |
| [PGSuper-Skew-System](PGSuper-Skew-System.md) | 21 s, export same as baseline, IDS same as baseline, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-System.usBridge-ids.txt) | 50 of 392 | 318 | 8 | 16 | 0 | 0 |
| [PGSuper-Skew-Bearings-System](PGSuper-Skew-Bearings-System.md) | 36 s, export same as baseline, IDS same as baseline, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-Bearings-System.usBridge-ids.txt) | 78 of 392 | 310 | 0 | 4 | 0 | 0 |
| [PGSuper-Skew-AlongPier-System](PGSuper-Skew-AlongPier-System.md) | 32 s, export same as baseline, IDS same as baseline, usBridge IDS: 8 of 599 specifications fail (see PGSuper-Skew-AlongPier-System.usBridge-ids.txt) | 79 of 392 | 309 | 0 | 4 | 0 | 0 |
| [PGSuper-Superelevation-System](PGSuper-Superelevation-System.md) | 29 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-Superelevation-System.usBridge-ids.txt) | 319 of 473 | 141 | 13 | 0 | 0 | 0 |
| [PGSuper-Superelevation-NoProfile-System](PGSuper-Superelevation-NoProfile-System.md) | 44 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-Superelevation-NoProfile-System.usBridge-ids.txt), removed IfcAlignmentVertical before the import | 319 of 473 | 141 | 13 | 0 | 0 | 0 |
| [PGSuper-ExteriorGirders-System](PGSuper-ExteriorGirders-System.md) | 39 s, export same as baseline, IDS same as baseline, usBridge IDS: 9 of 599 specifications fail (see PGSuper-ExteriorGirders-System.usBridge-ids.txt) | 306 of 473 | 154 | 13 | 0 | 0 | 0 |
| [PennDOT](PennDOT.md) | 23 s, no expected values: 72 of 108 values set by the import | | | | | | |
| [PennDOT-StandardTable](PennDOT-StandardTable.md) | 20 s, log check: 2 of 2 found, no expected values: 62 of 108 values set by the import | | | | | | |
| [Iowa](Iowa.md) | 64 s, no expected values: 230 of 274 values set by the import | | | | | | |
| [Iowa-StandardTable](Iowa-StandardTable.md) | 29 s, log check: 3 of 3 found, no expected values: 167 of 274 values set by the import | | | | | | |
| Standard table <-> IDS | OK: same: 326 declarations (material property sets not compared) | | | | | | |
| Standard table vs usBridge workbook | DIFFERS: 4 CLASS, 3 REQUIRED, 38 NO-VALUE, 13 VALUE, 54 TYPE, 1 NAME, 17 ENUM (20 PSET, 23 RESTRICTED, 36 OPTIONAL, 27 EXTRA, 1 NO-CLASS) (see Standard-vs-usBridge-workbook.txt) | | | | | | |
