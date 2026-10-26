# IFC import validation

- Run: 2026-10-26 08:23
- Configuration: Regression:Regression

| Run | Result | Match | Default match | Mismatch | Not imported | Not comparable | Missing |
|---|---|---|---|---|---|---|---|
| [PGSuper-Display](PGSuper-Display.md) | 8 s, export same as baseline, IDS same as baseline | 449 of 616 | 154 | 13 | 0 | 0 | 0 |
| [PGSuper-System](PGSuper-System.md) | 8 s, export same as baseline, IDS same as baseline | 449 of 616 | 154 | 13 | 0 | 0 | 0 |
| [PGSuper-WSDOT-System](PGSuper-WSDOT-System.md) | 5 s, export same as baseline, IDS same as baseline, configuration "WSDOT (Local)":"WSDOT (Local)" | 462 of 616 | 154 | 0 | 0 | 0 | 0 |
| [PGSuper-LongitudinalRebar-System](PGSuper-LongitudinalRebar-System.md) | 8 s, export same as baseline, IDS same as baseline | 449 of 616 | 154 | 13 | 0 | 0 | 0 |
| [PGSuper-Skew-System](PGSuper-Skew-System.md) | 9 s, export same as baseline, IDS same as baseline | 273 of 612 | 319 | 8 | 12 | 0 | 0 |
| [PGSuper-Skew-Bearings-System](PGSuper-Skew-Bearings-System.md) | 9 s, export same as baseline, IDS same as baseline | 301 of 612 | 311 | 0 | 0 | 0 | 0 |
| [PGSuper-Skew-AlongPier-System](PGSuper-Skew-AlongPier-System.md) | 10 s, export same as baseline, IDS same as baseline | 302 of 612 | 310 | 0 | 0 | 0 | 0 |
| [PGSuper-Superelevation-System](PGSuper-Superelevation-System.md) | 9 s, export same as baseline, IDS same as baseline | 462 of 616 | 141 | 13 | 0 | 0 | 0 |
| [PGSuper-Superelevation-NoProfile-System](PGSuper-Superelevation-NoProfile-System.md) | 8 s, export same as baseline, IDS same as baseline, removed IfcAlignmentVertical before the import | 462 of 616 | 141 | 13 | 0 | 0 | 0 |
| [PGSuper-Superelevation-Skew-System](PGSuper-Superelevation-Skew-System.md) | 9 s, export same as baseline, IDS same as baseline | 465 of 616 | 138 | 13 | 0 | 0 | 0 |
| [PGSuper-ExteriorGirders-System](PGSuper-ExteriorGirders-System.md) | 14 s, export same as baseline, IDS same as baseline | 449 of 616 | 154 | 13 | 0 | 0 | 0 |
| [PennDOT](PennDOT.md) | 6 s, no expected values: 78 of 108 values set by the import | | | | | | |
| [PennDOT-StandardTable](PennDOT-StandardTable.md) | 5 s, log check: 2 of 2 found, no expected values: 62 of 108 values set by the import | | | | | | |
| [Iowa](Iowa.md) | 27 s, no expected values: 233 of 274 values set by the import | | | | | | |
| [Iowa-StandardTable](Iowa-StandardTable.md) | 23 s, log check: 3 of 3 found, no expected values: 170 of 274 values set by the import | | | | | | |
| Standard table <-> IDS | OK: same: 340 declarations (material property sets not compared) | | | | | | |
