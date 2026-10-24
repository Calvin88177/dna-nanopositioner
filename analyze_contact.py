"""
analyze_contact.py <name> <h> <traj> [<skip_frames>]

Reads a contact trajectory and reports, in nm:
  d      axial distance between a reference slab of the stage (bp 4-11 above its face) and
         a reference slab of the post (bp 7-14 below its face). Its change between the
         piston and flush structures is the step; its variance gives the stiffness of
         everything between the two slabs (contact + piston tip + a few bp of bundle).
  tip    axial gap between the last base pair of the piston and the first base pair of the
         stage's centre helix (0.34 nm when coaxially stacked)
  face   axial gap between post face and stage face, averaged over the eight outer helices
  bp     whether the base pairs next to the contact are intact
"""
import sys
import os
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "tacoxDNA", "src"))
from oxtraj import frames, n_particles, block_sem, autocorr_time, SU_LEN_NM, KT_25C_PNNM
from setup_contact import load_map, ids_in, pad_of, FACE, STAGE_LEN

name, h, traj = sys.argv[1], int(sys.argv[2]), sys.argv[3]
skip = int(sys.argv[4]) if len(sys.argv) > 4 else 0
kT = float(sys.argv[5]) if len(sys.argv) > 5 else KT_25C_PNNM

m = load_map()
vs = json.load(open(name + ".json"))["vstrands"]
rc = {v["num"]: (v["row"], v["col"]) for v in vs}
centre = [n for n, x in rc.items() if x == (1, 1)][0]
outer = [n for n, x in rc.items() if x not in pad_of(name)]
stage_lo = FACE + h

g_stage = ids_in(m, stage_lo + 4, stage_lo + 12)
g_post = ids_in(m, FACE - 15, FACE - 7)


def pair(vh, vb):
    s = m._scaf.get((vh, vb))
    t = m._stap.get((vh, vb))
    if not s or not t or not s[1] or not t[1]:
        return None
    return s[1][0], t[1][0]


tip_post = pair(centre, FACE + h - 1)
tip_stage = pair(centre, stage_lo)
face_post = [pair(o, FACE - 1) for o in outer]
face_stage = [pair(o, stage_lo) for o in outer]
watch = {"piston tip bp": pair(centre, FACE + h - 1), "piston tip-1 bp": pair(centre, FACE + h - 2),
         "piston tip-2 bp": pair(centre, FACE + h - 3), "stage centre first bp": pair(centre, stage_lo),
         "stage centre second bp": pair(centre, stage_lo + 1)}

N = n_particles(name + ".json.top")
T, D, TIP, FACEGAP, LAT = [], [], [], [], []
BP = {k: [] for k in watch}
for t, a in frames(traj, N):
    pos, a1 = a[:, :3], a[:, 3:6]
    zdir = np.sign(pos[g_stage, 2].mean() - pos[g_post, 2].mean())
    T.append(t)
    D.append(zdir * (pos[g_stage, 2].mean() - pos[g_post, 2].mean()))
    zs = pos[list(tip_stage)].mean(0)
    zp = pos[list(tip_post)].mean(0)
    TIP.append(zdir * (zs[2] - zp[2]))
    LAT.append(np.linalg.norm(zs[:2] - zp[:2]))
    FACEGAP.append(np.mean([zdir * (pos[list(s)].mean(0)[2] - pos[list(p)].mean(0)[2])
                            for s, p in zip(face_stage, face_post)]))
    for k, pr in watch.items():
        i, j = pr
        # distance between the two base sites (0.4 along a1 from the centre of mass)
        BP[k].append(np.linalg.norm((pos[i] + 0.4 * a1[i]) - (pos[j] + 0.4 * a1[j])))

T = np.array(T)
sl = slice(skip, None)
print(f"{name}: {len(T)} frames, using {len(T[sl])} after skipping {skip}")
for lab, x in (("d (slab to slab)", D), ("tip gap", TIP), ("face gap", FACEGAP), ("tip lateral offset", LAT)):
    x = np.array(x)[sl] * SU_LEN_NM
    tau = autocorr_time(x)
    print(f"  {lab:<22s} mean {x.mean():7.3f} nm   sd {x.std(ddof=1):6.3f} nm   sem {block_sem(x):6.3f}   tau {tau:5.1f} frames")
d = np.array(D)[sl] * SU_LEN_NM
tip = np.array(TIP)[sl] * SU_LEN_NM
print(f"  stiffness from var(d):   {kT / d.var(ddof=1) / 1000:6.2f} N/m      from var(tip gap): {kT / tip.var(ddof=1) / 1000:6.2f} N/m")
for k, v in BP.items():
    v = np.array(v)[sl]
    print(f"  {k:<24s} base-base distance median {np.median(v):.2f} su, fraction open (>0.7 su): {np.mean(v > 0.7):.3f}")
np.save(os.path.splitext(traj)[0] + "_series.npy", np.array([T, np.array(D) * SU_LEN_NM, np.array(TIP) * SU_LEN_NM,
                                                               np.array(FACEGAP) * SU_LEN_NM]))
