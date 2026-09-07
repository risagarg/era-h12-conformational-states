#!/usr/bin/env python3
"""Distribution of H12 states across all human ERα LBD structures in the PDB.
Joint plot: 2D landscape + marginal histograms; plus a state-count bar."""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from collections import Counter

rows = [r for r in csv.DictReader(open('results/tables/sweep_results.csv'))
        if r['x'] not in ('', 'None') and r['y'] not in ('', 'None')]
data = [(float(r['x']), float(r['y']), r['state'], r['pdb']) for r in rows]
data = [d for d in data if d[0] <= 20 and d[1] <= 20]   # drop 2 artifacts

COL = {'agonist':'tab:green','antagonist':'firebrick','INTERMEDIATE':'darkorange','near-corner':'0.55'}
counts = Counter(d[2] for d in data)
n = len(data)

fig = plt.figure(figsize=(10, 9))
gs = GridSpec(2, 2, width_ratios=[4,1], height_ratios=[1,4], hspace=0.04, wspace=0.04)
ax = fig.add_subplot(gs[1,0])
axx = fig.add_subplot(gs[0,0], sharex=ax)   # top marginal (H12->3ERT)
axy = fig.add_subplot(gs[1,1], sharey=ax)   # right marginal (H12->1GWR)

# main scatter
for d in data:
    big = d[2] == 'INTERMEDIATE'
    ax.scatter(d[0], d[1], s=110 if big else 24, color=COL.get(d[2],'gray'),
               alpha=0.95 if big else 0.45, edgecolor='white' if big else 'none',
               linewidth=0.7, zorder=5 if big else 2)
for d in data:
    if d[2] == 'INTERMEDIATE':
        ax.annotate(d[3], (d[0], d[1]), fontsize=8, fontweight='bold', xytext=(6,3), textcoords='offset points')
ax.scatter(0,15.1, marker='*', s=300, color='black', zorder=6)
ax.scatter(15.1,0, marker='*', s=300, color='black', zorder=6)
ax.annotate('3ERT\nantagonist', (0,15.1), fontsize=8, xytext=(8,-22), textcoords='offset points')
ax.annotate('1GWR\nagonist', (15.1,0), fontsize=8, xytext=(-20,14), textcoords='offset points')
ax.set_xlabel('H12 backbone RMSD to 3ERT (antagonist)  →')
ax.set_ylabel('H12 backbone RMSD to 1GWR (agonist)  →')
ax.grid(alpha=0.25); ax.set_xlim(-1,17); ax.set_ylim(-1,17)

# marginal histograms
bins = [i*0.5 for i in range(0,35)]
axx.hist([d[0] for d in data], bins=bins, color='steelblue', edgecolor='white')
axy.hist([d[1] for d in data], bins=bins, orientation='horizontal', color='steelblue', edgecolor='white')
axx.tick_params(labelbottom=False); axy.tick_params(labelleft=False)
axx.set_ylabel('count'); axy.set_xlabel('count')
axx.set_title(f'H12 conformational states across {n} human ERα LBD structures in the PDB', fontsize=12, pad=10)
# shade the intermediate valley
for a, horiz in ((axx, False), (axy, True)):
    if horiz: a.axhspan(3, 12, color='orange', alpha=0.10)
    else:     a.axvspan(3, 12, color='orange', alpha=0.10)

# legend / counts box
order = ['agonist','antagonist','INTERMEDIATE','near-corner']
txt = "  ".join(f"{k}: {counts.get(k,0)}" for k in order if counts.get(k))
fig.text(0.5, 0.015, f"{txt}    (intermediate = {counts.get('INTERMEDIATE',0)}/{n} = "
         f"{100*counts.get('INTERMEDIATE',0)/n:.1f}%)", ha='center', fontsize=10)
fig.savefig('results/figures/h12_states_distribution.png', dpi=150, bbox_inches='tight')
print("state counts:", dict(counts), " n=", n)
print("wrote results/figures/h12_states_distribution.png")
