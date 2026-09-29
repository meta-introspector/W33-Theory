"""Pass 11136 scan (needs the Pass 11108 dumps and the Pass 11119 diag.py helpers): lowest neutral / charged tachyon levels on the diagonal T1 = T2 = B + i y for 36621 and 40521; frozen as data/w33_pass11136_domain_samples.json"""
import sys, json, re
exec(open('diag.py').read().split("def run(a):")[0])
from multiprocessing import Pool
MODELS = ['A8SM_20260942_36621_TF', 'A8SM_20260971_40521_TF']
PTS = [(0.0, 1.2), (0.0, 1.0), (0.25, 1.2), (0.5, 1.2), (0.5, 1.0), (0.5, 0.9), (0.5, 0.8), (0.4, 0.9), (0.3, 0.9),
       (0.3, 0.8), (0.45, 0.7), (0.2, 0.8), (0.35, 1.0), (0.6, 1.2), (0.8, 1.2), (1.0, 1.4), (1.0, 1.0), (0.75, 0.9),
       (1.5, 1.2), (1.5, 0.9), (0.5, 0.6), (0.1, 0.6)]
def job(a):
    lab, idx, rows, x, y = a
    US, TH = parsed(D2); us, th = US[idx], TH[idx]
    return lab, x, y, classes(us, th, rows, x + 1j * y, x + 1j * y)
if __name__ == '__main__':
    resc = json.load(open(R + r'\w33_pass11108_rescan_levels.json'))
    idx2 = {lab: int(i) for i, lab in (re.match(r'MODEL (\d+) (\S+)', l).groups() for l in open(D2[0]) if l.startswith('MODEL'))}
    jobs = [(lab, idx2[lab], resc[lab]['rows'], x, y) for lab in MODELS for x, y in PTS]
    with Pool(2) as p:
        res = p.map(job, jobs, chunksize=2)
    out = {}
    for lab, x, y, c in res:
        out.setdefault(lab, []).append(dict(B=x, y=y, **{k: float(v) for k, v in c.items()}))
        print(lab[14:], x, y, {k: round(float(v), 6) for k, v in c.items()}, flush=True)
    json.dump(out, open('domains.json', 'w'), indent=1)
