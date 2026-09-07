#!/usr/bin/env python3
"""
Tool-agnostic H12-state scorer for ERα LBD predictions.

Every prediction (AF2/ColabFold, Chai-1, Boltz-2; templated or cofolded) is
scored identically so all runs land as points on ONE landscape:

  axis X = H12 backbone RMSD to 3ERT  (antagonist / displaced H12)
  axis Y = H12 backbone RMSD to 1GWR  (agonist   / closed H12)

Both are computed AFTER superposing on the LBD core (305-530, which excludes
H12 at 537-547) so the RMSD measures *where H12 sits relative to the body of
the domain*, not global placement. A model in the antagonist corner has
low X / high Y; an agonist model has high X / low Y.

Auxiliary metrics (same for every run):
  - E523-K548 distance  (CA-CA, and sidechain carboxylate-Nz min)  -> state proxy
  - H12 %helix           (phi/psi-based)                            -> validity check

Usage:
  python score_h12.py MODEL.cif [MODEL2.pdb ...]
  python score_h12.py --validate          # score the two references themselves
"""
import sys, json, math
import gemmi

ANTAGONIST_REF = ('references/3ert.cif', 'A')   # 4-OHT bound, displaced H12
AGONIST_REF    = ('references/1gwr.cif', 'A')   # E2 bound, closed H12

CORE = (305, 530)          # LBD core for superposition (H12 @ 537-547 is outside)
H12  = (537, 547)          # the state-defining helix
BB   = ('N', 'CA', 'C', 'O')


def protein_chain(model, want=None):
    """Return the (largest) protein chain, or a named one."""
    best, best_n = None, -1
    for ch in model:
        poly = [r for r in ch if _is_aa(r)]
        if want is not None:
            if ch.name == want:
                return ch
            continue
        if len(poly) > best_n:
            best, best_n = ch, len(poly)
    return best


def _is_aa(r):
    info = gemmi.find_tabulated_residue(r.name)
    return bool(info and info.is_amino_acid())


def ca_map(chain, lo, hi):
    """{resnum: Position} for CA atoms in [lo,hi]."""
    out = {}
    for r in chain:
        if _is_aa(r) and lo <= r.seqid.num <= hi:
            a = r.find_atom('CA', '*')
            if a is not None:
                out[r.seqid.num] = a.pos
    return out


def bb_map(chain, lo, hi):
    """{(resnum,atom): Position} for backbone atoms in [lo,hi]."""
    out = {}
    for r in chain:
        if _is_aa(r) and lo <= r.seqid.num <= hi:
            for a in BB:
                at = r.find_atom(a, '*')
                if at is not None:
                    out[(r.seqid.num, a)] = at.pos
    return out


def superpose_core(model_chain, ref_chain):
    """Superpose model core CAs onto ref core CAs. Returns (Transform, core_rmsd, n)."""
    mm = ca_map(model_chain, *CORE)
    rm = ca_map(ref_chain, *CORE)
    common = sorted(set(mm) & set(rm))
    if len(common) < 4:
        return None, float('nan'), len(common)
    fixed  = [rm[n] for n in common]   # reference stays put
    moving = [mm[n] for n in common]   # model is moved onto it
    sup = gemmi.superpose_positions(fixed, moving)
    return sup.transform, sup.rmsd, len(common)


def h12_rmsd_after(model_chain, ref_chain, transform):
    """RMSD of model H12 backbone (after `transform`) to ref H12 backbone."""
    mh = bb_map(model_chain, *H12)
    rh = bb_map(ref_chain, *H12)
    keys = sorted(set(mh) & set(rh))
    if not keys:
        return float('nan'), 0
    s = 0.0
    for k in keys:
        p = transform.apply(mh[k]) if transform is not None else mh[k]
        q = rh[k]
        s += (p.x-q.x)**2 + (p.y-q.y)**2 + (p.z-q.z)**2
    return math.sqrt(s/len(keys)), len(keys)


def _dist(a, b):
    return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)


# NOTE: position 548 is ARG in canonical P03372 (the "K548" label is a misnomer;
# the residue is arginine). So this is the E523-R548 carboxylate<->guanidinium pair.
ACID_ATOMS = {'GLU': ('OE1', 'OE2'), 'ASP': ('OD1', 'OD2')}
BASE_ATOMS = {'ARG': ('NH1', 'NH2', 'NE'), 'LYS': ('NZ',)}


def e523_548(chain):
    """E523<->R548 distances: CA-CA, and min(acidic O -> basic N) salt-bridge dist.
    Robust to the actual residue type at each position."""
    a523 = b548 = None
    for r in chain:
        if r.seqid.num == 523:
            a523 = r
        if r.seqid.num == 548:
            b548 = r
    out = {'ca_ca': None, 'salt_bridge': None, 'res523': None, 'res548': None}
    if a523 and b548:
        out['res523'], out['res548'] = a523.name, b548.name
        aca, bca = a523.find_atom('CA', '*'), b548.find_atom('CA', '*')
        if aca and bca:
            out['ca_ca'] = round(_dist(aca.pos, bca.pos), 2)
        acids = [a523.find_atom(n, '*') for n in ACID_ATOMS.get(a523.name, ())]
        bases = [b548.find_atom(n, '*') for n in BASE_ATOMS.get(b548.name, ())]
        ds = [_dist(o.pos, n.pos) for o in acids if o for n in bases if n]
        if ds:
            out['salt_bridge'] = round(min(ds), 2)
    return out


def h12_pct_helix(chain):
    """Fraction of H12 residues with backbone phi/psi in the alpha-helical basin."""
    res = [r for r in chain if _is_aa(r) and H12[0] <= r.seqid.num <= H12[1]]
    res_all = [r for r in chain if _is_aa(r)]
    by_num = {r.seqid.num: r for r in res_all}
    def atom(r, name):
        a = r.find_atom(name, '*'); return a.pos if a else None
    helical = 0; counted = 0
    for r in res:
        prev = by_num.get(r.seqid.num - 1); nxt = by_num.get(r.seqid.num + 1)
        if not (prev and nxt): continue
        C0 = atom(prev, 'C'); N = atom(r,'N'); CA = atom(r,'CA'); C = atom(r,'C'); N1 = atom(nxt,'N')
        if None in (C0, N, CA, C, N1): continue
        phi = math.degrees(gemmi.calculate_dihedral(C0, N, CA, C))
        psi = math.degrees(gemmi.calculate_dihedral(N, CA, C, N1))
        counted += 1
        if -120 <= phi <= -30 and -80 <= psi <= -5:
            helical += 1
    return round(helical / counted, 2) if counted else None


def score(model_path, model_chain_name=None, offset=0):
    """offset is added to every model residue number before scoring, to map a
    tool's 1..N output onto canonical P03372 numbering (construct 155-595 -> +154).
    References are always in canonical numbering and untouched."""
    st = gemmi.read_structure(model_path)
    mc = protein_chain(st[0], model_chain_name)
    if offset:
        for r in mc:
            r.seqid.num += offset
    row = {'model': model_path, 'chain': mc.name, 'offset': offset}
    for label, (ref_path, ref_cn) in (('3ert_antag', ANTAGONIST_REF),
                                      ('1gwr_agon', AGONIST_REF)):
        rst = gemmi.read_structure(ref_path)
        rc = protein_chain(rst[0], ref_cn)
        T, core_rmsd, n_core = superpose_core(mc, rc)
        h12r, n_h12 = h12_rmsd_after(mc, rc, T)
        row[f'h12rmsd_{label}'] = None if math.isnan(h12r) else round(h12r, 2)
        row[f'core_rmsd_{label}'] = None if math.isnan(core_rmsd) else round(core_rmsd, 2)
        row[f'n_core_{label}'] = n_core
        row[f'n_h12_{label}'] = n_h12
    row['e523_548'] = e523_548(mc)
    row['h12_pct_helix'] = h12_pct_helix(mc)
    return row


if __name__ == '__main__':
    args = sys.argv[1:]
    offset = 0
    if '--offset' in args:
        i = args.index('--offset'); offset = int(args[i + 1]); del args[i:i + 2]
    if not args or args[0] == '--validate':
        targets = [('references/3ert.cif', 'A'), ('references/1gwr.cif', 'A')]
        print("=== VALIDATION: scoring the two reference structures ===")
        print("(expect each ref ~0 to itself, large to the other; opposite corners)\n")
        rows = [score(p, c) for p, c in targets]
    else:
        rows = [score(p, offset=offset) for p in args]
    for r in rows:
        print(json.dumps(r))
    print("\n# axes: h12rmsd_3ert_antag (X)  vs  h12rmsd_1gwr_agon (Y)")
