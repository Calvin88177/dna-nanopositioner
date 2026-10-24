"""
setup_bend.py - three-point bending of a free square-lattice blade in oxDNA.

Why three-point bending: a blade of a parallelogram flexure is a fixed-guided beam, whose
deflected shape is antisymmetric about mid-span. Each half is a cantilever of length L/2.
A simply supported beam of the same length loaded at mid-span has the same two halves, with
the "clamp" supplied by symmetry (zero slope and zero inter-layer slip at the centre), and
the free ends correspond to the blade's inflection point. So with central load P and
mid-span sag s, the two-blade parallelogram stiffness for blades of length L is

        k_parallelogram = (P/2) / s          (Euler-Bernoulli: 24 B / L^3)

Supports: every nucleotide in one base-pair plane near each end is held in the thin (y)
direction only by a harmonic trap (oxDNA "twist" force with zero rate and a mask).
Load: a constant force along y shared by the nucleotides of the two central planes.

usage: setup_bend.py <prefix> <length_bp> <support_lo> <support_hi> <load_lo> <load_hi> <P_pN> <tag> <conf>
"""
import pickle
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tacoxDNA", "src"))

SU_FORCE_PN = 48.63      # 1 oxDNA force unit in pN
SU_LEN_NM = 0.8518       # 1 oxDNA length unit in nm


def load_map(path="virt2nuc"):
    with open(path, "rb") as f:
        vh_vb2nuc, pattern = pickle.load(f)
    return vh_vb2nuc


def plane_ids(m, vbases):
    ids = []
    for (vh, vb), (strand, nucs) in list(m._scaf.items()) + list(m._stap.items()):
        if vb in vbases:
            ids.extend(nucs)
    return sorted(set(ids))


def read_conf(path):
    with open(path) as f:
        lines = f.readlines()
    box = [float(x) for x in lines[1].split("=")[1].split()]
    data = np.array([[float(x) for x in l.split()[:9]] for l in lines[3:] if l.strip()])
    return box, data


def write_forces(path, conf, supports, load, P_su, stiff=2.0, hold_xz=0.02):
    """supports: list of id-lists; load: id-list; P_su: total central force along +y (oxDNA units)."""
    box, data = read_conf(conf)
    pos = data[:, :3]
    out = []
    for k, ids in enumerate(supports):
        for i in ids:
            x, y, z = pos[i]
            # knife-edge support: restrain y only
            out.append("{\ntype = twist\nparticle = %d\nstiff = %g\nrate = 0\nbase = 0\n"
                       "pos0 = %.6f, %.6f, %.6f\ncenter = %.6f, %.6f, %.6f\naxis = 0, 0, 1\nmask = 0, 1, 0\n}\n"
                       % (i, stiff, x, y, z, x, y, z))
            # very weak in-plane hold so the blade does not drift or spin in the x-z plane;
            # the first support holds x and z, the second only x, so the span is free to shorten
            mask = "1, 0, 1" if k == 0 else "1, 0, 0"
            out.append("{\ntype = twist\nparticle = %d\nstiff = %g\nrate = 0\nbase = 0\n"
                       "pos0 = %.6f, %.6f, %.6f\ncenter = %.6f, %.6f, %.6f\naxis = 0, 0, 1\nmask = %s\n}\n"
                       % (i, hold_xz, x, y, z, x, y, z, mask))
    if abs(P_su) > 0:
        f_each = P_su / len(load)
        ids = ",".join(str(i) for i in load)
        out.append("{\ntype = string\nparticle = %s\nF0 = %.8f\nrate = 0\ndir = 0, %d, 0\n}\n"
                   % (ids, abs(f_each), 1 if f_each > 0 else -1))
    with open(path, "w") as f:
        f.write("\n".join(out))


INPUT = """backend = CPU
backend_precision = double
sim_type = MD
interaction_type = DNA2
salt_concentration = 1.0
use_average_seq = 1
T = {T}
dt = 0.005
verlet_skin = 0.05
thermostat = brownian
newtonian_steps = 103
diff_coeff = 2.5
steps = {steps}
topology = {top}
conf_file = {conf}
trajectory_file = {tag}_traj_full.dat
lastconf_file = {tag}_last.dat
energy_file = {tag}_energy.dat
print_conf_interval = {steps}
print_energy_every = 10000
time_scale = linear
refresh_vel = 1
restart_step_counter = 1
no_stdout_energy = 1
external_forces = 1
external_forces_file = {tag}_forces.txt
log_file = {tag}.log
data_output_1 = {{
    name = {tag}_traj.dat
    print_every = {every}
    col_1 = {{
        type = configuration
        print_momenta = false
    }}
}}
"""

if __name__ == "__main__":
    prefix, length, s_lo, s_hi, l_lo, l_hi, P_pN, tag, conf = sys.argv[1:10]
    steps = sys.argv[10] if len(sys.argv) > 10 else "1000000"
    every = sys.argv[11] if len(sys.argv) > 11 else "1000"
    T = sys.argv[12] if len(sys.argv) > 12 else "25C"
    m = load_map()
    A = plane_ids(m, {int(s_lo)})
    B = plane_ids(m, {int(s_hi)})
    M = plane_ids(m, {int(l_lo), int(l_hi)})
    print(f"support planes: {len(A)} and {len(B)} nucleotides; load planes: {len(M)} nucleotides")
    P_su = float(P_pN) / SU_FORCE_PN
    write_forces(f"{tag}_forces.txt", conf, [A, B], M, P_su)
    with open(f"{tag}_input", "w") as f:
        f.write(INPUT.format(T=T, steps=steps, top=prefix + ".top", conf=conf, tag=tag, every=every))
    np.save(f"{tag}_groups.npy", np.array([A, B, M], dtype=object), allow_pickle=True)
    print(f"wrote {tag}_input, {tag}_forces.txt  (P = {P_pN} pN = {P_su:.4f} su)")
