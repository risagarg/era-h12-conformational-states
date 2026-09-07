#!/usr/bin/env python3
"""
Verify candidate ERα LBD PDBs by structure content (not memory) and score their
H12 on our landscape. An ID is accepted only if its protein chain matches
canonical P03372 numbering (so it's human ERα LBD); ligands are read from the file.
"""
import os, urllib.request, math
import gemmi, score_h12 as S

CANON = "".join(l.strip() for l in open('inputs/P03372.fasta') if not l.startswith('>'))
LIG_NAMES = {  # CCD -> readable / mode hint (only for annotation)
    'EST': 'estradiol (agonist)', 'DES': 'diethylstilbestrol (agonist)',
    'OHT': '4-hydroxytamoxifen (SERM)', 'RAL': 'raloxifene (SERM)',
    'GW5': 'GW5638 (SERD/carboxylate)',
}
# broad candidate set across modes; verification filters wrong/non-ER IDs
CANDIDATES = ["1GWR","1ERE","1A52","1QKU","3ERD","1GWQ","3ERT","1ERR","2OUZ","1XPC",
              "1XP1","1XP6","1SJ0","2BJ4","4XI3","5W9C","1R5K","5T1Z","6CHW","5ACC",
              "2IOK","5UFX","6B0F","2P15"]

os.makedirs('refs', exist_ok=True)

def fetch(pid):
    # files pre-downloaded via curl into references/candidates/<PID>.cif
    for p in (f'references/candidates/{pid}.cif', f'references/candidates/{pid.lower()}.cif'):
        if os.path.exists(p) and os.path.getsize(p) > 3000:
            return p
    return None

def is_aa(r):
    i = gemmi.find_tabulated_residue(r.name); return bool(i and i.is_amino_acid())

print(f"{'PDB':<6}{'ER?':<5}{'span':<12}{'H12->3ERT':<11}{'H12->1GWR':<11}{'ligands'}")
print('-'*78)
rows = []
for pid in CANDIDATES:
    path = fetch(pid)
    if not path:
        print(f"{pid:<6}  (not fetched / 404)"); continue
    try:
        st = gemmi.read_structure(path)
    except Exception:
        print(f"{pid:<6}  (parse error)"); continue
    prot = None
    for ch in st[0]:
        poly = [r for r in ch if is_aa(r)]
        if len(poly) > 150:
            nums = {r.seqid.num: r for r in poly}
            if 540 in nums and 353 in nums:
                prot = ch; break
    if prot is None:
        print(f"{pid:<6} no   (no ER-LBD-like chain)"); continue
    poly = [r for r in prot if is_aa(r)]
    mm = sum(1 for r in poly if 1 <= r.seqid.num <= len(CANON)
             and gemmi.find_tabulated_residue(r.name).one_letter_code.upper() != CANON[r.seqid.num-1])
    is_er = mm <= 3
    nums = [r.seqid.num for r in poly]
    ligs = sorted({r.name for ch in st[0] for r in ch
                   if not is_aa(r) and r.name not in ('HOH','CCS','SO4','GOL','EDO','CL','NA','ACT','PEG')})
    h3 = h1 = None
    if is_er:
        for label, (rp, rc) in (('3ert', S.ANTAGONIST_REF), ('1gwr', S.AGONIST_REF)):
            rst = gemmi.read_structure(rp); refc = S.protein_chain(rst[0], rc)
            T, _, _ = S.superpose_core(prot, refc)
            h, n = S.h12_rmsd_after(prot, refc, T)
            if label == '3ert': h3 = None if math.isnan(h) else round(h, 1)
            else: h1 = None if math.isnan(h) else round(h, 1)
    ligstr = ", ".join(f"{l}={LIG_NAMES.get(l,'?')}" for l in ligs) if ligs else "(none)"
    print(f"{pid:<6}{'yes' if is_er else 'no':<5}{f'{min(nums)}-{max(nums)}':<12}{str(h3):<11}{str(h1):<11}{ligstr}")
    if is_er:
        rows.append((pid, h3, h1, ligs))

print(f"\nverified ER LBD structures: {len(rows)}")
