"""
analyze_bend.py <dir> <prefix> <tag> [<skip_frames>]

Sag of a three-point-bending run: y of the loaded mid-span planes minus the mean y of the
two support planes, per frame. Reports mean, spread, block standard error and correlation time.
"""
import sys
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from oxtraj import frames, n_particles, block_sem, autocorr_time, SU_LEN_NM, KT_25C_PNNM


def sag_series(d, prefix, tag):
    A, B, M = np.load(os.path.join(d, f"{tag}_groups.npy"), allow_pickle=True)
    A, B, M = list(A), list(B), list(M)
    N = n_particles(os.path.join(d, prefix + ".top"))
    T, S, SPAN = [], [], []
    for t, a in frames(os.path.join(d, f"{tag}_traj.dat"), N):
        pos = a[:, :3]
        ya, yb, ym = pos[A, 1].mean(), pos[B, 1].mean(), pos[M, 1].mean()
        T.append(t)
        S.append((ym - 0.5 * (ya + yb)) * SU_LEN_NM)
        SPAN.append(abs(pos[B, 2].mean() - pos[A, 2].mean()) * SU_LEN_NM)
    return np.array(T), np.array(S), np.array(SPAN)


if __name__ == "__main__":
    d, prefix, tag = sys.argv[1:4]
    skip = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    T, S, SPAN = sag_series(d, prefix, tag)
    s = S[skip:]
    print(f"{tag}: {len(S)} frames, using {len(s)} after skipping {skip}; last step {T[-1]:.0f}")
    print(f"  span between supports   {SPAN[skip:].mean():6.2f} nm")
    print(f"  sag    mean {s.mean():7.3f} nm   sd {s.std(ddof=1):6.3f} nm   sem {block_sem(s):6.3f} nm   tau {autocorr_time(s):5.1f} frames")
    print(f"  stiffness of the sag coordinate from its variance: {KT_25C_PNNM / s.var(ddof=1) / 1000:6.3f} N/m"
          f"   -> two-blade equivalent (half of it): {KT_25C_PNNM / s.var(ddof=1) / 2000:6.3f} N/m")
    np.save(os.path.join(d, f"{tag}_sag.npy"), np.array([T, S]))
