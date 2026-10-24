"""analyze_duplex.py <tag> [<skip>] - height and base-pair state of the two-duplex contact."""
import sys
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "tacoxDNA", "src"))
from oxtraj import frames, block_sem, autocorr_time, SU_LEN_NM, KT_25C_PNNM
from setup_duplex import load_map, ids_in

tag = sys.argv[1]
skip = int(sys.argv[2]) if len(sys.argv) > 2 else 50
m = load_map()
anchor, handle = ids_in(m, 0, 8), ids_in(m, 24, 32)


def pair(vb):
    return m._scaf[(0, vb)][1][0], m._stap[(0, vb)][1][0]


watch = {vb: pair(vb) for vb in (12, 13, 14, 15, 16, 17, 18, 19)}
T, D, GAP = [], [], []
OPEN = {vb: [] for vb in watch}
for t, a in frames(f"{tag}_traj.dat", 64):
    pos, a1 = a[:, :3], a[:, 3:6]
    zdir = np.sign(pos[handle, 2].mean() - pos[anchor, 2].mean())
    T.append(t)
    D.append(zdir * (pos[handle, 2].mean() - pos[anchor, 2].mean()) * SU_LEN_NM)
    GAP.append(zdir * (pos[list(watch[16])].mean(0)[2] - pos[list(watch[15])].mean(0)[2]) * SU_LEN_NM)
    for vb, (i, j) in watch.items():
        OPEN[vb].append(np.linalg.norm((pos[i] + 0.4 * a1[i]) - (pos[j] + 0.4 * a1[j])) > 0.7)
D, GAP = np.array(D), np.array(GAP)
d, g = D[skip:], GAP[skip:]
print(f"{tag}: {len(D)} frames (last step {T[-1]:.0f}), using {len(d)}")
print(f"  slab-to-slab height  mean {d.mean():6.3f} nm  sd {d.std(ddof=1):5.3f}  sem {block_sem(d):5.3f}  "
      f"stiffness {KT_25C_PNNM / d.var(ddof=1) / 1000:5.2f} N/m   tau {autocorr_time(d):4.1f} frames")
print(f"  contact gap          mean {g.mean():6.3f} nm  sd {g.std(ddof=1):5.3f}")
print("  fraction of frames with the base pair open:  " +
      "  ".join(f"bp{vb}: {np.mean(OPEN[vb][skip:]):.2f}" for vb in watch))
q = np.array_split(d, 4)
print("  height by quarter of the run: " + "  ".join(f"{x.mean():.3f}" for x in q))
np.save(f"{tag}_series.npy", np.array([T, D, GAP]))
