# IFC import validation

- Run: 2026-10-15 11:16
- Configuration: Regression:Regression

| Run | Result | Match | Default match | Mismatch | Not imported | Not comparable | Missing |
|---|---|---|---|---|---|---|---|
| [PGSuper-Display](PGSuper-Display.md) | 3 s, export differs from baseline: 1482 lines (see PGSuper-Display.export-diff.txt), IDS differs from baseline: 780 lines (see PGSuper-Display.ids-diff.txt) | 124 of 239 | 101 | 0 | 14 | 0 | 0 |
| [PGSuper-System](PGSuper-System.md) | 3 s, export differs from baseline: 1482 lines (see PGSuper-System.export-diff.txt), IDS differs from baseline: 780 lines (see PGSuper-System.ids-diff.txt) | 124 of 239 | 101 | 0 | 14 | 0 | 0 |
| [PGSuper-Skew-System](PGSuper-Skew-System.md) | 6 s, export differs from baseline: 2280 lines (see PGSuper-Skew-System.export-diff.txt), IDS differs from baseline: 1200 lines (see PGSuper-Skew-System.ids-diff.txt) | 30 of 294 | 220 | 8 | 36 | 0 | 0 |
| [PGSuper-Skew-Bearings-System](PGSuper-Skew-Bearings-System.md) | 5 s, export differs from baseline: 2280 lines (see PGSuper-Skew-Bearings-System.export-diff.txt), IDS differs from baseline: 1200 lines (see PGSuper-Skew-Bearings-System.ids-diff.txt) | 58 of 294 | 212 | 0 | 24 | 0 | 0 |
| [PGSuper-Skew-AlongPier-System](PGSuper-Skew-AlongPier-System.md) | 5 s, export differs from baseline: 2280 lines (see PGSuper-Skew-AlongPier-System.export-diff.txt), IDS differs from baseline: 1200 lines (see PGSuper-Skew-AlongPier-System.ids-diff.txt) | 59 of 294 | 211 | 0 | 24 | 0 | 0 |
| [PennDOT](PennDOT.md) | 2 s, no expected values: 59 of 102 values set by the import | | | | | | |
| [PennDOT-StandardTable](PennDOT-StandardTable.md) | 2 s, log check: 2 of 2 found, no expected values: 45 of 102 values set by the import | | | | | | |
| [Iowa](Iowa.md) | 8 s, no expected values: 178 of 260 values set by the import | | | | | | |
| [Iowa-StandardTable](Iowa-StandardTable.md) | 7 s, log check: 3 of 3 found, no expected values: 118 of 260 values set by the import | | | | | | |
| Standard table <-> IDS | OK: same: 272 declarations (material property sets not compared) | | | | | | |
