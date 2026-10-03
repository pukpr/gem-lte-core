import json, sys
from concurrent.futures import ThreadPoolExecutor
from three import make3, score, B
JOBS = [("s3_s25_B2.5_W4", "s25", 2.5, 4.0), ("s3_s25_B2.5_W4.2", "s25", 2.5, 4.2), ("s3_amo_B0.2319_W4.2", "amo", 0.231889340633804, 4.2)]
import single_k
def go(j):
    name, s, kB, kW = j; d = make3(name, s, kB, kW)
    x = score(d, search=300); print(json.dumps(x), flush=True); return x
with ThreadPoolExecutor(3) as ex: recs = list(ex.map(go, JOBS))
json.dump(recs, open(B / "three_search.json", "w"), indent=1); print("DONE")
