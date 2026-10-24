"""blade_profile.py <dir> <json> <traj> - inter-helix spacing statistics of a blade over a trajectory."""
import sys
import os
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "tacoxDNA", "src"))
from oxtraj import frames, n_particles, SU_LEN_NM
import pickle

d, js, traj = sys.argv[1], sys.argv[2], sys.argv[3]
lo, hi = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (34, 94)
os.chdir(d)
with open("virt2nuc", "rb") as f:
    m, _ = pickle.load(f)
vs = json.load(open(js))["vstrands"]
rc2num = {(v["row"], v["col"]): v["num"] for v in vs}
rows = max(r for r, c in rc2num) + 1
cols = max(c for r, c in rc2num) + 1
N = n_particles(js + ".top")


def axis(pos, vh, vb):
    s = m._scaf[(vh, vb)][1][0]
    t = m._stap[(vh, vb)][1][0]
    return 0.5 * (pos[s] + pos[t])


print("step   in-layer spacing (nm): mean / max    inter-layer spacing (nm): mean / max    [span bp %d-%d]" % (lo, hi - 1))
for t, a in frames(traj, N):
    pos = a[:, :3]
    inl, intl = [], []
    for vb in range(lo, hi, 2):
        ax = {rc: axis(pos, n, vb) for rc, n in rc2num.items()}
        for r in range(rows):
            for c in range(cols - 1):
                inl.append(np.linalg.norm(ax[(r, c)] - ax[(r, c + 1)]))
        for c in range(cols):
            for r in range(rows - 1):
                intl.append(np.linalg.norm(ax[(r, c)] - ax[(r + 1, c)]))
    inl = np.array(inl) * SU_LEN_NM
    intl = np.array(intl) * SU_LEN_NM
    print("%7d   %5.2f / %5.2f                        %5.2f / %5.2f" % (t, inl.mean(), inl.max(), intl.mean(), intl.max()))
