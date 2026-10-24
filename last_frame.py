"""last_frame.py <top> <traj> <out> - write the last complete frame of a trajectory as a full configuration."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oxtraj import n_particles

top, traj, out = sys.argv[1:4]
n = n_particles(top)
with open(traj) as f:
    lines = f.readlines()
starts = [i for i, l in enumerate(lines) if l.startswith("t =")]
for s in reversed(starts):
    block = lines[s:s + 3 + n]
    if len(block) == 3 + n and all(len(l.split()) >= 9 for l in block[3:]):
        with open(out, "w") as g:
            g.writelines(block[:3])
            for l in block[3:]:
                p = l.split()[:9]
                g.write(" ".join(p) + " 0 0 0 0 0 0\n")
        print("wrote", out, "from frame", block[0].strip())
        break
