#!/usr/bin/env python3
"""
Score H12 of EVERY human ERα LBD structure on the landscape and find all
intermediate (off-both-corner) H12 states. Refs loaded once for speed.
"""
import glob, os, math, csv, json
import gemmi, score_h12 as S
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CANON = "".join(l.strip() for l in open('inputs/P03372.fasta') if not l.startswith('>'))
SKIP = {'HOH','SO4','GOL','EDO','CL','NA','ACT','PEG','MG','ZN','BME','TRS','IOD','DMS',
        'CCS','ACY','PO4','MES','EPE','FLC','CIT','K','CA','MN','NO3','FMT','DMSO','BME'}

def is_aa(r):
    i = gemmi.find_tabulated_residue(r.name); return bool(i and i.is_amino_acid())

# load references once
ra = S.protein_chain(gemmi.read_structure('references/3ert.cif')[0], 'A')
rg = S.protein_chain(gemmi.read_structure('references/1gwr.cif')[0], 'A')

def classify(x, y):
    if x is None or y is None: return 'no-H12'
    if x > 25 or y > 25:       return 'artifact'
    if y <= 2.5 and x >= 7:    return 'agonist'
    if x <= 2.5 and y >= 7:    return 'antagonist'
    if x >= 3 and y >= 3:      return 'INTERMEDIATE'
    return 'near-corner'

rows = []
for path in sorted(glob.glob('references/pdb_sweep/*.cif')):
    pid = os.path.basename(path)[:-4].upper()
    if os.path.getsize(path) < 3000: continue
    try:
        st = gemmi.read_structure(path)
    except Exception:
        continue
    prot = None
    for ch in st[0]:
        poly = [r for r in ch if is_aa(r)]
        if len(poly) > 150 and any(r.seqid.num == 353 for r in ch) and any(r.seqid.num == 540 for r in ch):
            prot = ch; break
    if prot is None:
        continue  # not an LBD-with-H12 structure (e.g. DBD, or H12 unmodeled)
    poly = [r for r in prot if is_aa(r)]
    mm = sum(1 for r in poly if 1 <= r.seqid.num <= len(CANON)
             and gemmi.find_tabulated_residue(r.name).one_letter_code.upper() != CANON[r.seqid.num-1])
    if mm > 15:
        continue  # not ERα
    try:
        T1, _, _ = S.superpose_core(prot, ra); x, nx = S.h12_rmsd_after(prot, ra, T1)
        T2, _, _ = S.superpose_core(prot, rg); y, ny = S.h12_rmsd_after(prot, rg, T2)
    except Exception:
        continue
    x = None if (x is None or math.isnan(x)) else round(x, 1)
    y = None if (y is None or math.isnan(y)) else round(y, 1)
    ligs = sorted({r.name for c in st[0] for r in c if not is_aa(r) and r.name not in SKIP and len(list(r)) >= 8})
    rows.append(dict(pdb=pid, x=x, y=y, state=classify(x, y), muts=mm, ligs=";".join(ligs)))

from collections import Counter
dist = Counter(r['state'] for r in rows)
print(f"scored {len(rows)} ERα LBD structures (with modeled H12)")
print("state distribution:", dict(dist))
inter = [r for r in rows if r['state'] == 'INTERMEDIATE']
inter.sort(key=lambda r: r['x'])
print(f"\n=== {len(inter)} INTERMEDIATE H12 structures (off both corners) ===")
print(f"  {'PDB':<6}{'H12->3ERT':<11}{'H12->1GWR':<11}{'ligands'}")
for r in inter:
    print(f"  {r['pdb']:<6}{str(r['x']):<11}{str(r['y']):<11}{r['ligs']}")

with open('results/tables/sweep_results.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['pdb','x','y','state','muts','ligs']); w.writeheader(); w.writerows(rows)

# plot all, highlight intermediates
fig, ax = plt.subplots(figsize=(8.5, 7.5))
col = {'agonist':'tab:green','antagonist':'firebrick','INTERMEDIATE':'orange',
       'near-corner':'0.6','artifact':'0.85'}
for r in rows:
    if r['x'] is None or r['y'] is None: continue
    big = r['state'] == 'INTERMEDIATE'
    ax.scatter(r['x'], r['y'], s=90 if big else 22, color=col.get(r['state'],'gray'),
               alpha=0.95 if big else 0.4, edgecolor='white' if big else 'none',
               linewidth=0.6, zorder=5 if big else 2)
    if big:
        ax.annotate(r['pdb'], (r['x'], r['y']), fontsize=7, xytext=(4,3), textcoords='offset points')
ax.scatter(0,15.1, marker='*', s=300, color='black'); ax.scatter(15.1,0, marker='*', s=300, color='black')
ax.annotate('3ERT (antag)', (0,15.1), fontsize=8, xytext=(6,-2), textcoords='offset points')
ax.annotate('1GWR (agon)', (15.1,0), fontsize=8, xytext=(6,4), textcoords='offset points')
ax.set_xlabel('H12 RMSD to 3ERT (antagonist)'); ax.set_ylabel('H12 RMSD to 1GWR (agonist)')
ax.set_title(f'All human ERα LBD H12 states (n={len([r for r in rows if r["x"] is not None])})\n'
             f'orange = intermediate (n={len(inter)})')
ax.grid(alpha=0.25); fig.tight_layout(); fig.savefig('results/figures/sweep_landscape.png', dpi=150)
print("\nwrote results/tables/sweep_results.csv, results/figures/sweep_landscape.png")
