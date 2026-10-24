"""
setup_contact.py - a piston-on-stage contact under constant preload, in oxDNA.

Two coaxial 3x3 square-lattice bundles:
    post   helices occupy [0, FACE); the centre helix (the piston) runs to FACE + h
    stage  all nine helices occupy [FACE + h, FACE + h + 48)
so the stage sits end-on on the piston (h > 0) or on the post face (h = 0).

Guidance and load (standing in for the flexure, which is not simulated):
    post, lowest 24 bp     held in x, y and z   (the anchor)
    stage, highest 24 bp   held in x and y only, and pushed towards the post with a
                           constant total force (the preload)
Everything within about 8 nm of the contact is free.

Steps:
    setup_contact.py design <name> <h>                    write caDNAno JSON
    setup_contact.py inputs <name> <h> <tag> <conf> <F_pN> <trap_stiff> <steps> [<every>] [<T>]
"""
import sys
import os
import pickle
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tacoxDNA", "src"))
sys.path.insert(0, HERE)

FACE = 54                 # post face; with margin 6 the piston has its last crossover at 47|48, 8 bp below a 2-bp tip
STAGE_LEN = 48
HOLD = 24                 # base pairs held at each far end
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
topology = {name}.json.top
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


# which helices of the post stand proud: the centre helix alone, or a 2x2 pad that includes it
PADS = {"c4": [(0, 0), (0, 1), (1, 0), (1, 1)]}


def pad_of(name):
    return PADS.get(name, [(1, 1)])


def segments(h, name=""):
    hel = {}
    for r in range(3):
        for c in range(3):
            top = FACE + h if (r, c) in pad_of(name) else FACE
            hel[(r, c)] = [(0, top), (FACE + h, FACE + h + STAGE_LEN)]
    return hel


def load_map(path="virt2nuc"):
    with open(path, "rb") as f:
        m, pattern = pickle.load(f)
    return m


def ids_in(m, lo, hi, vhs=None):
    out = []
    for (vh, vb), (strand, nucs) in list(m._scaf.items()) + list(m._stap.items()):
        if lo <= vb < hi and (vhs is None or vh in vhs):
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
    mode, name, h = sys.argv[1], sys.argv[2], int(sys.argv[3])
    if mode == "design":
        from make_block import write_cadnano
        write_cadnano(segments(h, name), name + ".json", margin=6)
    else:
        tag, conf, F_pN, stiff, steps = sys.argv[4], sys.argv[5], float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
        every = sys.argv[9] if len(sys.argv) > 9 else "200"
        T = sys.argv[10] if len(sys.argv) > 10 else "25C"
        guide = sys.argv[11] if len(sys.argv) > 11 else "far"       # "far" or "ideal"
        free = int(sys.argv[12]) if len(sys.argv) > 12 else 6
        m = load_map()
        stage_lo, stage_hi = FACE + h, FACE + h + STAGE_LEN
        anchor = ids_in(m, 0, HOLD)
        handle = ids_in(m, stage_hi - HOLD, stage_hi)
        pos = read_conf(conf)[:, :3]
        zdir = np.sign(pos[handle, 2].mean() - pos[anchor, 2].mean())      # stage sits at +zdir
        out = [twist(i, stiff, pos[i], "1, 1, 1") for i in anchor]
        if guide == "far":
            out += [twist(i, stiff, pos[i], "1, 1, 0") for i in handle]
            n_lat = len(handle)
        else:
            # ideal guide: every nucleotide further than `free` bp from the contact is held in x and y
            vs = json.load(open(name + ".json"))["vstrands"]
            proud = [v["num"] for v in vs if (v["row"], v["col"]) in pad_of(name)]
            lat = []
            for (vh, vb), (strand, nucs) in list(m._scaf.items()) + list(m._stap.items()):
                if vb >= stage_lo:
                    held = vb >= stage_lo + free
                else:
                    top = FACE + h if vh in proud else FACE
                    held = HOLD <= vb < top - free
                if held:
                    lat.extend(nucs)
            lat = sorted(set(lat))
            out += [twist(i, stiff, pos[i], "1, 1, 0") for i in lat]
            n_lat = len(lat)
        if F_pN != 0:
            f_each = abs(F_pN) / SU_FORCE_PN / len(handle)
            sgn = -zdir if F_pN > 0 else zdir                              # positive = pressing on
            out.append("{\ntype = string\nparticle = %s\nF0 = %.8f\nrate = 0\ndir = 0, 0, %d\n}\n"
                       % (",".join(str(i) for i in handle), f_each, int(sgn)))
        with open(f"{tag}_forces.txt", "w") as f:
            f.write("\n".join(out))
        with open(f"{tag}_input", "w") as f:
            f.write(INPUT.format(T=T, steps=steps, name=name, conf=conf, tag=tag, every=every))
        print(f"{tag}: guide {guide}, anchor {len(anchor)} nt, laterally held {n_lat} nt, handle {len(handle)} nt, "
              f"stage towards z{'+' if zdir > 0 else '-'}, F = {F_pN} pN, trap stiffness {stiff}")
