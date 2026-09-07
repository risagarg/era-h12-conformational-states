#!/usr/bin/env python3
"""Build the ground-truth ERα H12 landscape from verified crystal structures.
Each structure scored on our two axes; grouped by ligand binding-mode class."""
import os, math, csv, glob
import gemmi, score_h12 as S
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CANON = "".join(l.strip() for l in open('inputs/P03372.fasta') if not l.startswith('>'))

# resolved ligand CCD -> (readable, mode). Mode = ligand class (H12 is measured, not assumed).
LIG = {
 'EST': ('estradiol','agonist'), 'DES': ('diethylstilbestrol','agonist'),
 'EZT': ('estradiol-vinyl deriv','agonist'), 'Q97': ('ethoxyphenyl-butene diphenol','agonist'),
 'ZTW': ('raloxifene CORE (no side chain)','agonist'),
 'OHT': ('4-hydroxytamoxifen','SERM'), 'RAL': ('raloxifene','SERM'),
 'C3D': ('lasofoxifene-type','SERM'), 'AIT': ('benzoxathiin SERM','SERM'),
 'AIH': ('benzoxathiin SERM','SERM'), 'AIU': ('benzoxathiin SERM','SERM'),
 'E4D': ('benzoxathiin SERM','SERM'), '29S': ('bazedoxifene','SERM'),
 'IOK': ('indole-acetamide SERM','SERM'), '86Y': ('benzopyran SERM','SERM'),
 'F3D': ('indazole-amide (SERD-like)','SERM'),
 'GW5': ('GW5638','SERD'), 'KE9': ('acrylic-acid SERD','SERD'), 'C6V': ('LSZ102','SERD'),
}
SKIP = {'HOH','CCS','SO4','GOL','EDO','CL','NA','ACT','PEG','AU','ACY','DMS','MG','ZN','BME','TRS','IOD'}
MODE_COLOR = {'agonist':'tab:green','SERM':'firebrick','SERD':'darkorange','?':'gray'}

def is_aa(r):
    i = gemmi.find_tabulated_residue(r.name); return bool(i and i.is_amino_acid())

rows = []
for path in sorted(glob.glob('references/candidates/*.cif')):
    pid = os.path.basename(path)[:-4].upper()
    st = gemmi.read_structure(path)
    prot = None
    for ch in st[0]:
        if len([r for r in ch if is_aa(r)]) > 150 and any(r.seqid.num == 353 for r in ch):
            prot = ch; break
    if prot is None:
        continue
    poly = [r for r in prot if is_aa(r)]
    muts = [r.seqid.num for r in poly if 1 <= r.seqid.num <= len(CANON)
            and gemmi.find_tabulated_residue(r.name).one_letter_code.upper() != CANON[r.seqid.num-1]]
    if len(muts) > 12:
        continue  # not ERα
    ligs = [r.name for c in st[0] for r in c if not is_aa(r) and r.name not in SKIP]
    lig = next((l for l in ligs if l in LIG), ligs[0] if ligs else None)
    name, mode = LIG.get(lig, (lig or 'apo', '?'))
    # score H12 on landscape
    sc = {}
    for label, (rp, rc) in (('3ert', S.ANTAGONIST_REF), ('1gwr', S.AGONIST_REF)):
        rst = gemmi.read_structure(rp); refc = S.protein_chain(rst[0], rc)
        T, _, _ = S.superpose_core(prot, refc)
        h, n = S.h12_rmsd_after(prot, refc, T)
        sc[label] = (None if math.isnan(h) else round(h, 1), n)
    x, nx = sc['3ert']; y, ny = sc['1gwr']
    rows.append(dict(pdb=pid, lig=lig, name=name, mode=mode, muts=len(muts),
                     n_h12=nx, x=x, y=y))

# table + csv
rows.sort(key=lambda r: (r['mode'], r['x'] if r['x'] is not None else 99))
print(f"{'PDB':<6}{'lig':<5}{'mode':<9}{'H12->3ERT':<11}{'H12->1GWR':<11}{'muts':<5}{'ligand'}")
print('-'*82)
for r in rows:
    print(f"{r['pdb']:<6}{str(r['lig']):<5}{r['mode']:<9}{str(r['x']):<11}{str(r['y']):<11}{r['muts']:<5}{r['name']}")
with open('results/tables/reference_map.csv','w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=['pdb','lig','name','mode','muts','n_h12','x','y']); w.writeheader(); w.writerows(rows)

# plot
fig, ax = plt.subplots(figsize=(8.5, 7.5))
ax.scatter(0, 15.1, marker='*', s=380, color='black', zorder=6); ax.annotate('3ERT', (0,15.1), fontsize=8, xytext=(6,4), textcoords='offset points')
ax.scatter(15.1, 0, marker='*', s=380, color='black', zorder=6); ax.annotate('1GWR', (15.1,0), fontsize=8, xytext=(6,4), textcoords='offset points')
seen=set()
for r in rows:
    if r['x'] is None or r['y'] is None: continue
    m=r['mode']; lbl=m if m not in seen else None; seen.add(m)
    ax.scatter(r['x'], r['y'], color=MODE_COLOR.get(m,'gray'), s=130, alpha=0.85, edgecolor='white', label=lbl, zorder=4)
    ax.annotate(f"{r['pdb']}", (r['x'], r['y']), fontsize=7, xytext=(5,3), textcoords='offset points')
ax.set_xlabel('H12 RMSD to 3ERT (antagonist)  ->  agonist-like')
ax.set_ylabel('H12 RMSD to 1GWR (agonist)  ->  antagonist-like')
ax.set_title('Ground-truth ERα H12 landscape (verified crystals, by ligand mode)')
ax.legend(); ax.grid(alpha=0.25); fig.tight_layout(); fig.savefig('results/figures/reference_map.png', dpi=150)
print(f"\n{len(rows)} structures mapped -> results/tables/reference_map.csv, results/figures/reference_map.png")
PYEOF_GUARD = None
