"""oxtraj.py - minimal reader for oxDNA text trajectories and a block-average helper."""
import numpy as np

SU_LEN_NM = 0.8518
SU_FORCE_PN = 48.63
KT_25C_PNNM = 4.1164          # kB * 298.15 K in pN nm


def frames(path, n):
    """yield (time, array[n, ncols]) for each frame of a trajectory with n particles."""
    with open(path) as f:
        while True:
            head = f.readline()
            if not head:
                return
            if not head.startswith("t"):
                continue
            t = float(head.split("=")[1])
            f.readline()
            f.readline()
            rows = [f.readline() for _ in range(n)]
            if len(rows) < n or not rows[-1].strip():
                return                                  # incomplete last frame
            try:
                arr = np.array([r.split()[:9] for r in rows], dtype=float)
            except ValueError:
                return
            if arr.shape != (n, 9):
                return
            yield t, arr


def n_particles(top):
    with open(top) as f:
        return int(f.readline().split()[0])


def block_sem(x, nblocks=5):
    """standard error of the mean from block averages (correlated samples)."""
    x = np.asarray(x, dtype=float)
    if len(x) < 2 * nblocks:
        return float("nan")
    b = np.array([c.mean() for c in np.array_split(x, nblocks)])
    return b.std(ddof=1) / np.sqrt(nblocks)


def autocorr_time(x):
    """integrated autocorrelation time in units of samples (simple windowed estimate)."""
    x = np.asarray(x, dtype=float) - np.mean(x)
    n = len(x)
    if n < 8 or x.var() == 0:
        return float("nan")
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    tau = 0.5
    for k in range(1, n):
        if acf[k] <= 0:
            break
        tau += acf[k]
    return tau
