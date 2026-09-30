# IFC import validation

- Run: 2026-09-29 20:43
- Configuration: Regression:Regression

| Run | Result | Match | Default match | Mismatch | Not imported | Not comparable | Missing |
|---|---|---|---|---|---|---|---|
| [PGSuper-Display](PGSuper-Display.md) | 28 s, export same as baseline, IDS same as baseline | 124 of 239 | 101 | 0 | 14 | 0 | 0 |
| [PGSuper-System](PGSuper-System.md) | 35 s, export same as baseline, IDS same as baseline | 124 of 239 | 101 | 0 | 14 | 0 | 0 |
| [PGSuper-Skew-System](PGSuper-Skew-System.md) | 40 s, export same as baseline, IDS same as baseline | 30 of 294 | 220 | 8 | 36 | 0 | 0 |
| [PGSuper-Skew-Bearings-System](PGSuper-Skew-Bearings-System.md) | 41 s, export same as baseline, IDS same as baseline | 58 of 294 | 212 | 0 | 24 | 0 | 0 |
| [PGSuper-Skew-AlongPier-System](PGSuper-Skew-AlongPier-System.md) | 38 s, export same as baseline, IDS same as baseline | 59 of 294 | 147 | 64 | 24 | 0 | 0 |
| [PennDOT](PennDOT.md) | 12 s, no expected values: 59 of 102 values set by the import | | | | | | |
| [PennDOT-StandardTable](PennDOT-StandardTable.md) | 11 s, log check: 2 of 2 found, no expected values: 45 of 102 values set by the import | | | | | | |
| [Iowa](Iowa.md) | 24 s, no expected values: 178 of 260 values set by the import | | | | | | |
| [Iowa-StandardTable](Iowa-StandardTable.md) | 22 s, log check: 3 of 3 found, no expected values: 118 of 260 values set by the import | | | | | | |
| Standard table <-> IDS | OK: same: 272 declarations (material property sets not compared) | | | | | | |
