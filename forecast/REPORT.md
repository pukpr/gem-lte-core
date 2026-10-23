# GEM-LTE forecast ledger

Updated 2026-10-16 00:05 UTC. Latest frozen projection per target. `ss_*` = 1 - MSE(model)/MSE(baseline); positive beats the baseline. post_data = after the committed .dat ends (never fitted); post_issue = after the ledger timestamp (strictly ex-ante).

`obs cal r` = correlation of the live feed with the committed .dat over their overlap; LOW (< 0.9) means the live feed is not the same series and post-data scores for that target are not yet meaningful.

| target | issue | data end | obs cal r | r holdout | post-data n | r | SS persist | SS AR1 | SS clim | post-issue n | r | SS persist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| soi | 2f4a34bc3fe7 | 2026.5833 | 0.98 | 0.4385 | 0 |  |  |  |  | 0 |  |  |
| qbo30 | d5c551229093 | 2026.6667 | 1.00 | 0.8431 | 0 |  |  |  |  | 0 |  |  |
| qbo50 | af5577f5fd56 | 2024.2500 | 1.00 | 0.6503 | 28 | 0.4664 | 0.7218 | 0.4629 | 0.1065 | 0 |  |  |
| qbo70 | 0dd8b8f15fbc | 2026.6667 | 1.00 | 0.5029 | 0 |  |  |  |  | 0 |  |  |
| pdo | 3cda1196a7dc | 2022.0833 | 0.65 LOW | 0.3328 | 54 | 0.3413 | -7.8931 | -1.0668 | -0.1735 | 0 |  |  |
| nino4 | 2669440ecfdf | 2024.1667 | 0.97 | 0.8234 | 29 | -0.0899 | 0.6030 | 0.3571 | -0.0588 | 0 |  |  |
| nino34 | 2f946a43a8d0 | 2024.2500 | 0.99 | 0.4846 | 28 | -0.5795 | 0.1178 | -0.1049 | -0.5718 | 0 |  |  |
| baltic | 30a87d53f921 | 2025.5000 | 0.41 LOW | 0.0109 | 4 | -0.4704 | -220.2666 | -236.3985 | -1066.1393 | 0 |  |  |
| lod | 2026-10 | 2026.7036 |  |  | 0 |  |  |  |  | 0 |  |  |
| chandler_x | 2026-10 | 2026.7036 |  |  | 0 |  |  |  |  | 0 |  |  |
| chandler_y | 2026-10 | 2026.7036 |  |  | 0 |  |  |  |  | 0 |  |  |

## Group summary (median over members)

| group | members | median post-data r | median SS AR1 | share SS AR1 > 0 |
|---|---|---|---|---|
| qbo50 | 1 | 0.466 | 0.463 | 1.00 |
| pdo | 1 | 0.341 | -1.067 | 0.00 |
| nino4 | 1 | -0.090 | 0.357 | 1.00 |
| nino34 | 1 | -0.580 | -0.105 | 0.00 |
| baltic | 1 | -0.470 | -236.399 | 0.00 |

## Rolling-origin hindcasts (refit engines)

Each origin uses only data available at that origin; scored over the full horizon that followed.

| target | origins | median r | median SS persist | median SS AR1 | median SS clim | share SS AR1 > 0 |
|---|---|---|---|---|---|---|
| lod | 41 | 0.751 | 0.202 | 0.487 | 0.548 | 0.90 |
| chandler_x | 41 | 0.921 | 0.692 | 0.604 | 0.286 | 0.88 |
| chandler_y | 41 | 0.926 | 0.650 | 0.516 | -2.139 | 0.83 |
