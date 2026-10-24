"""
setup_duplex.py - two coaxial blunt-ended duplexes pressed end to end with a constant force.

The smallest model of a single-helix stop: a 16-bp "post" duplex (bp 0-15) and a 16-bp
"stage" duplex (bp 16-31) on the same helix axis, so the two blunt ends stack.
    post,  bp 0-7      held in x, y and z
    stage, bp 24-31    held in x and y, pushed towards the post with total force F
    bp 8-23            held in x and y except the last `free` bp on each side of the contact
                       (with 8 bp free on each side the two ends simply slide past each other)
Question: does the contact hold its height under load, or do the terminal base pairs
get squeezed out so that the stage sinks by one base-pair step at a time?

    setup_duplex.py design
    setup_duplex.py inputs <tag> <conf> <F_pN> <steps> <every> <T>
"""
import sys
import os
import pickle
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tacoxDNA", "src"))
sys.path.insert(0, HERE)
SU_FORCE_PN = 48.63

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
topology = dup.json.top
conf_file = {conf}
trajectory_file = {tag}_traj_full.dat
lastconf_file = {tag}_last.dat
energy_file = {tag}_energy.dat
print_conf_interval = {steps}
print_energy_every = 100000
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


def load_map(path="virt2nuc"):
    with open(path, "rb") as f:
        m, pattern = pickle.load(f)
    return m


def ids_in(m, lo, hi):
    out = []
    for (vh, vb), (strand, nucs) in list(m._scaf.items()) + list(m._stap.items()):
        if lo <= vb < hi:
            out.extend(nucs)
    return sorted(set(out))


def read_conf(path):
    with open(path) as f:
        lines = f.readlines()
    return np.array([[float(x) for x in l.split()[:9]] for l in lines[3:] if l.strip()])


def twist(i, stiff, p, mask):
    return ("{\ntype = twist\nparticle = %d\nstiff = %g\nrate = 0\nbase = 0\n"
            "pos0 = %.6f, %.6f, %.6f\ncenter = %.6f, %.6f, %.6f\naxis = 0, 0, 1\nmask = %s\n}\n"
            % (i, stiff, p[0], p[1], p[2], p[0], p[1], p[2], mask))


if __name__ == "__main__":
    if sys.argv[1] == "design":
        from make_block import write_cadnano
        write_cadnano({(0, 0): [(0, 16), (16, 32)]}, "dup.json")
    else:
        tag, conf, F_pN, steps, every, T = sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7]
        m = load_map()
        anchor = ids_in(m, 0, 8)
        handle = ids_in(m, 24, 32)
        pos = read_conf(conf)[:, :3]
        zdir = np.sign(pos[handle, 2].mean() - pos[anchor, 2].mean())
        free = int(sys.argv[8]) if len(sys.argv) > 8 else 3
        guide = ids_in(m, 8, 16 - free) + ids_in(m, 16 + free, 24)     # sideways guidance up to `free` bp from the contact
        out = [twist(i, 2.0, pos[i], "1, 1, 1") for i in anchor]
        out += [twist(i, 2.0, pos[i], "1, 1, 0") for i in handle + guide]
        if F_pN != 0:
            out.append("{\ntype = string\nparticle = %s\nF0 = %.8f\nrate = 0\ndir = 0, 0, %d\n}\n"
                       % (",".join(str(i) for i in handle), abs(F_pN) / SU_FORCE_PN / len(handle),
                          int(-zdir if F_pN > 0 else zdir)))
        with open(f"{tag}_forces.txt", "w") as f:
            f.write("\n".join(out))
        with open(f"{tag}_input", "w") as f:
            f.write(INPUT.format(T=T, steps=steps, conf=conf, tag=tag, every=every))
        print(f"{tag}: F = {F_pN} pN at {T}, {steps} steps")
