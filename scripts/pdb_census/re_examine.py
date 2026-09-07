#!/usr/bin/env python3
"""Rigorous re-check: score EVERY chain, flag bad superpositions, show the full
distribution + threshold sensitivity, and a metric-independent H12-centroid view."""
import glob, os, math, csv
import gemmi, score_h12 as S

CANON = "".join(l.strip() for l in open('inputs/P03372.fasta') if not l.startswith('>'))
def is_aa(r):
    i = gemmi.find_tabulated_residue(r.name); return bool(i and i.is_amino_acid())

ra = S.protein_chain(gemmi.read_structure('references/3ert.cif')[0], 'A')
rg = S.protein_chain(gemmi.read_structure('references/1gwr.cif')[0], 'A')
# 3ERT H12 centroid in its own frame (common frame = 3ERT core)
def h12_centroid(chain):
    ps = [a.pos for r in chain if is_aa(r) and 537 <= r.seqid.num <= 547
          for a in r if a.name == 'CA']
    if not ps: return None
    return gemmi.Position(sum(p.x for p in ps)/len(ps), sum(p.y for p in ps)/len(ps), sum(p.z for p in ps)/len(ps))
ra_h12c = h12_centroid(ra); rg_T,_,_ = S.superpose_core(rg, ra)  # 1gwr core onto 3ert frame
# 1gwr H12 centroid mapped into 3ert frame
rg_h12c_raw = h12_centroid(rg); rg_h12c = rg_T.apply(rg_h12c_raw)

ids = [l.strip() for l in open('references/p03372_ids.txt')]
chains = []   # (pid, chain, x, y, core_rmsd, h12cen_in_3ert_frame)
for pid in ids:
    p = f'references/pdb_sweep/{pid}.cif'
    if not os.path.exists(p) or os.path.getsize(p) < 3000: continue
    try: st = gemmi.read_structure(p)
    except: continue
    for ch in st[0]:
        poly = [r for r in ch if is_aa(r)]
        nums = {r.seqid.num for r in poly}
        if len(poly) < 150 or 353 not in nums: continue
        if sum(1 for n in range(537,548) if n in nums) < 8: continue
        mm = sum(1 for r in poly if 1<=r.seqid.num<=len(CANON)
                 and gemmi.find_tabulated_residue(r.name).one_letter_code.upper()!=CANON[r.seqid.num-1])
        if mm > 15: continue
        try:
            T1,c1,_ = S.superpose_core(ch, ra); x,_ = S.h12_rmsd_after(ch, ra, T1)
            T2,c2,_ = S.superpose_core(ch, rg); y,_ = S.h12_rmsd_after(ch, rg, T2)
            cen = T1.apply(h12_centroid(ch))   # this chain's H12 centroid in 3ert frame
        except: continue
        if x is None or y is None or math.isnan(x) or math.isnan(y): continue
        chains.append((pid, ch.name, round(x,1), round(y,1), round(c1,1), cen))

print(f"scored {len(chains)} CHAINS (all chains, not just first) from {len(set(c[0] for c in chains))} structures")
bad = [c for c in chains if c[4] > 3.0]
print(f"  chains with poor core superposition (core RMSD>3, x/y unreliable): {len(bad)}")
good = [c for c in chains if c[4] <= 3.0 and c[2] <= 20 and c[3] <= 20]

# fine distribution along H12->3ERT axis
print("\n=== distribution of H12->3ERT RMSD (fine bins) ===")
import collections
b = collections.Counter(min(int(c[2]), 16) for c in good)
for k in range(0,17):
    bar = '#' * (b.get(k,0)//3)
    print(f"  {k:>2}-{k+1:<2} A: {b.get(k,0):>3} {bar}")

# threshold sensitivity for 'intermediate' (off BOTH corners)
print("\n=== how many INTERMEDIATE chains at different cutoffs (>= t from BOTH refs) ===")
for t in (2.0, 2.5, 3.0, 4.0):
    inter = [c for c in good if c[2] >= t and c[3] >= t]
    print(f"  cutoff {t} A: {len(inter)} chains  {sorted(set((c[0],c[1]) for c in inter))[:20]}")

# metric-independent: spread of H12 centroid in common (3ERT) frame
import statistics as stt
print("\n=== metric-independent: H12 centroid distance to the two reference centroids (common frame) ===")
def d(a,b): return math.sqrt((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2)
near_ant=near_ag=between=0
for c in good:
    da=d(c[5],ra_h12c); dg=d(c[5],rg_h12c)
    if da<4: near_ant+=1
    elif dg<4: near_ag+=1
    else: between+=1
print(f"  H12 centroid within 4A of 3ERT-H12: {near_ant}; within 4A of 1GWR-H12: {near_ag}; BETWEEN/elsewhere: {between}")

with open('results/tables/reexamine_chains.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['pdb','chain','x_to3ert','y_to1gwr','core_rmsd'])
    for c in chains: w.writerow(c[:5])
print("\nwrote results/tables/reexamine_chains.csv")
