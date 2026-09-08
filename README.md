# ERα Helix-12 Conformational-State Analysis

**Methods for measuring the conformational state of a nuclear-receptor switch, plus a reusable state-scoring approach.**

Helix 12 (H12) is the ~11-residue helix that acts as the on/off switch of the estrogen-receptor-α (ERα) ligand-binding domain. in one position it completes the coactivator-binding surface, allowing transcription of ER mediated genes, and in another it blocks the transcription. This repository is a **methods & skills extract** from a larger study that contains the reusable pieces; how to census H12 states across every deposited ERα structure, and how to score the H12 state of any structure.

This repo is methodology only, and it shows how the conformational-state analysis is done, does not contain any of the MD trajectories.

---

Files : 

1. PDB-wide H12 state census (`scripts/pdb_census/`)
Scan human ERα ligand-binding-domain structures in the PDB, superpose each onto a common rigid core (helices excluding H12), and classify its H12 as **agonist / antagonist / intermediate** by backbone RMSD to reference states.

![H12 state census](figures/sweep_landscape.png)
![H12 state distribution](figures/h12_states_distribution.png)

- `sweep_intermediates.py` - score every entry with a two-axis map (RMSD-to-agonist × RMSD-to-antagonist).
- `re_examine.py` — per-chain re-scoring (handles asymmetric dimers).
- `reference_map.py` — place the reference structures on the same map.
- `make_h12_histogram.py` — state distribution with marginals.
- `verify_refs.py` — sequence-verify candidate structures and read their bound ligand.

2. H12 conformational-state scoring (`scripts/h12_scoring/`)
Score the H12 state of any structure - Kabsch-superpose the ligand-binding-domain core (305–530 Cα, excluding H12), then measure H12 backbone RMSD to the agonist (1GWR) and antagonist (3ERT) references, plus geometric proxies — the E523–D545 opening distance and H12 % helicity. 

<p align="center">
  <img src="figures/fig_af2_groove.png" width="47%"/>
  <img src="figures/fig_af2_groove_antagonist_white.png" width="47%"/>
</p>

*The two states the scorer separates; agonist (left — coactivator peptide bound, H12 clear of the groove) vs antagonist (right — H12 folded into the coactivator groove).*

- `score_h12.py` — the core scorer (superposition + RMSD + geometry proxies).
- `analyze_cofold.py` — apply the scorer to predicted/cofolded structures with pose QC.



## Repo layout
```
scripts/
  pdb_census/    census H12 states across all public ERα PDB structures
  h12_scoring/   score H12 state of any structure (crystal / predicted / MD frame)
references/      reference crystal structures (1GWR agonist, 3ERT antagonist, 1R5K intermediate, 1ERR)
figures/         result figures on public data
```


- **Environment:** [`requirements.txt`](requirements.txt) (`gemmi`, `matplotlib`; PyMOL separately for renders).
- **Data:** `references/` holds the key structures. The full census set is every human-ERα LBD entry in the **public PDB** — downloadable by the accessions the census scripts iterate over; point the scripts at a local mirror to reproduce.
- **Conventions:** ERα is numbered by UniProt **P03372**; H12 = residues **537–547**; the superposition core = **305–530** (excludes H12); E523–R548/D545 opening is a clean state proxy (~10 Å agonist vs ~28–37 Å antagonist).

