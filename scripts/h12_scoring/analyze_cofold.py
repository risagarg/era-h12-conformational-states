#!/usr/bin/env python3
"""
Score every cofold pose with the validated two-axis landscape + 3 QC gates.

Per pose:
  axes        : H12 RMSD to 3ERT (antag) and to 1GWR (agon), LBD-core superposed
  GATE 1 pose : phenol O -> Glu353/Arg394 anchor (phenol ligands);
                carboxylate O -> Asp351/H11 (GW5638). Ligand-in-pocket check.
  GATE 2 helix: H12 %helix (well-formed helix, not coil)
  GATE 3 conf : Boltz confidence_score / iptm / ptm

Protein output is renumbered +154 (construct 155-595 -> canonical) before scoring.
"""
import sys, json, glob, math, os
import gemmi
import score_h12 as S

OFFSET = 154
# pocket anchor residues (canonical numbering)
GLU353 = ('GLU', 353, ('OE1', 'OE2'))
ARG394 = ('ARG', 394, ('NH1', 'NH2', 'NE'))
ASP351 = ('ASP', 351, ('OD1', 'OD2'))
H11 = (490, 530)         # H11 region for GW5638 carboxylate locus
POCKET = (353, 394, 521, 524, 351)   # for ligand-in-pocket centroid


def _is_aa(r):
    info = gemmi.find_tabulated_residue(r.name)
    return bool(info and info.is_amino_acid())


def load_renumbered(path):
    st = gemmi.read_structure(path)
    prot = S.protein_chain(st[0])
    for r in prot:
        if _is_aa(r):
            r.seqid.num += OFFSET
    # ligand chain = any chain that isn't the protein chain, with non-AA residues
    lig_atoms = []
    for ch in st[0]:
        if ch.name == prot.name:
            continue
        for r in ch:
            if not _is_aa(r):
                for a in r:
                    lig_atoms.append((a.element.name, a.pos))
    return st, prot, lig_atoms


def _dist(a, b):
    return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)


def res_atoms(chain, resname, resnum, names):
    for r in chain:
        if r.seqid.num == resnum:
            return [r.find_atom(n, '*') for n in names if r.find_atom(n, '*')]
    return []


def centroid(positions):
    n = len(positions)
    return gemmi.Position(sum(p.x for p in positions)/n,
                          sum(p.y for p in positions)/n,
                          sum(p.z for p in positions)/n)


def pose_qc(prot, lig_atoms, is_carboxylate_anchor):
    """Return dict of pose-sanity distances. Ligand atoms = [(element, pos)]."""
    lig_O = [pos for el, pos in lig_atoms if el == 'O']
    lig_all = [pos for el, pos in lig_atoms]
    out = {'n_lig_atoms': len(lig_all), 'n_lig_O': len(lig_O)}
    if not lig_all:
        out['ligand_present'] = False
        return out
    out['ligand_present'] = True

    # ligand-in-pocket: centroid distance ligand <-> pocket CAs
    pocket_ca = []
    for r in prot:
        if r.seqid.num in POCKET:
            a = r.find_atom('CA', '*')
            if a:
                pocket_ca.append(a.pos)
    if pocket_ca and lig_all:
        out['lig_to_pocket_centroid'] = round(_dist(centroid(lig_all), centroid(pocket_ca)), 2)

    def min_O_to(resspec):
        rn, rnum, names = resspec
        atoms = res_atoms(prot, rn, rnum, names)
        if not atoms or not lig_O:
            return None
        return round(min(_dist(o, a.pos) for o in lig_O for a in atoms), 2)

    if is_carboxylate_anchor:
        # GW5638: carboxylate O near Asp351 / H11
        out['ligO_to_Asp351'] = min_O_to(ASP351)
        # nearest ligand O to any H11 backbone CA
        h11_ca = [r.find_atom('CA', '*').pos for r in prot
                  if H11[0] <= r.seqid.num <= H11[1] and r.find_atom('CA', '*')]
        if h11_ca and lig_O:
            out['ligO_to_H11_min'] = round(min(_dist(o, c) for o in lig_O for c in h11_ca), 2)
        out['anchor_ok'] = (out.get('ligO_to_Asp351') is not None and out['ligO_to_Asp351'] < 6.0) \
                           or (out.get('ligO_to_H11_min') is not None and out['ligO_to_H11_min'] < 6.0)
    else:
        # phenol ligands: phenol O H-bonding Glu353/Arg394
        out['ligO_to_Glu353'] = min_O_to(GLU353)
        out['ligO_to_Arg394'] = min_O_to(ARG394)
        best = [d for d in (out['ligO_to_Glu353'], out['ligO_to_Arg394']) if d is not None]
        out['anchor_ok'] = bool(best) and min(best) < 3.5   # H-bond distance
    return out


def score_pose(path, carboxylate=False):
    st, prot, lig = load_renumbered(path)
    row = {'pose': os.path.basename(path)}
    for label, (ref_path, ref_cn) in (('3ert_antag', S.ANTAGONIST_REF),
                                      ('1gwr_agon', S.AGONIST_REF)):
        rst = gemmi.read_structure(ref_path)
        rc = S.protein_chain(rst[0], ref_cn)
        T, core_rmsd, n_core = S.superpose_core(prot, rc)
        h12r, n_h12 = S.h12_rmsd_after(prot, rc, T)
        row[f'h12_{label}'] = None if math.isnan(h12r) else round(h12r, 2)
        row[f'core_{label}'] = None if math.isnan(core_rmsd) else round(core_rmsd, 2)
    row['h12_pct_helix'] = S.h12_pct_helix(prot)
    row['qc'] = pose_qc(prot, lig, carboxylate)
    # confidence sidecar
    conf_path = path.replace('.cif', '.confidence.json')
    if os.path.exists(conf_path):
        c = json.load(open(conf_path))
        row['confidence'] = round(c.get('confidence_score', -1), 3)
        row['iptm'] = round(c.get('iptm', -1), 3) if isinstance(c.get('iptm'), (int, float)) else c.get('iptm')
    return row


if __name__ == '__main__':
    paths = sys.argv[1:] or sorted(glob.glob('runs/boltz/cofold_out/*.cif'))
    for p in paths:
        carbox = 'GW5638' in p
        r = score_pose(p, carboxylate=carbox)
        print(json.dumps(r))
