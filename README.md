# oxDNA component checks for the single-axis DNA nanopositioner

These scripts rebuild the brief's estimates and test two pieces of the device in oxDNA:
the guide blade and the stop contact. They are component models, not the device. The
full-device check still needs the caDNAno design.

Everything was run with oxDNA (CPU build, `DNA2` interaction, 1 M salt, average-sequence
parameters, 25 C, dt = 0.005, Brownian thermostat) and the tacoxDNA caDNAno converter.

## Files

| File | What it does |
| --- | --- |
| `flexure_estimates_check.py` | Independent rebuild of the brief's numbers (Part A) and the effects beam theory leaves out (Part B). Pure Python, no dependencies. |
| `make_block.py` | Writes a caDNAno square-lattice JSON for a block of parallel helices: one scaffold strand per helix, all staple crossovers at canonical positions, staples nicked into 32-48 nt pieces. |
| `setup_bend.py`, `analyze_bend.py` | Three-point bending of a 2x8 blade. Supports and load are oxDNA external forces. |
| `setup_contact.py`, `analyze_contact.py` | A 3x3 "post" with a protruding centre helix (the piston) pressed end-on against a 3x3 "stage" with a constant force. Two guidance modes: `far` (soft) and `ideal`. |
| `setup_duplex.py`, `analyze_duplex.py` | The smallest contact model: two coaxial duplexes pressed end to end, as a force series. |
| `oxtraj.py`, `last_frame.py`, `blade_profile.py` | Trajectory reader and helpers. |
| `RESULTS.md` | The numbers from the runs of 4-5 October 2026. |
| `designs/` | The caDNAno JSON files of the five test structures. |
| `run_*.sh`, `queue_*.sh` | The exact run sequences used. Paths are absolute; edit `OX` and `SIM` first. |

## Setup

- Build oxDNA (github.com/lorenzo-rovigatti/oxDNA) and set `OX` in the `.sh` files to the binary.
- Clone tacoxDNA (github.com/lorenzo-rovigatti/tacoxDNA) beside this folder, as `../tacoxDNA`.
  The Python scripts import its caDNAno index map from there.
- Python 3 with numpy.

The blade run, start to finish:

    python3 make_block.py 2 8 128 blade128.json
    python3 ../tacoxDNA/src/cadnano_oxDNA.py blade128.json sq -p     # writes .top, .oxdna and virt2nuc
    oxDNA relax_mc && oxDNA relax_md && oxDNA equil                  # inputs in inputs_blade128/
    python3 last_frame.py blade128.json.top equil_traj.dat start.dat
    python3 setup_bend.py blade128.json 128 34 93 63 64 70 Pp70 start.dat 120000 500
    oxDNA Pp70_input
    python3 analyze_bend.py . blade128.json Pp70 20

The contact and two-duplex runs are scripted whole in `run_contact_ideal.sh` and
`queue_duplex.sh`.

## Why three-point bending stands in for the parallelogram

A blade of a parallelogram flexure is a fixed-guided beam. Its deflected shape is
antisymmetric about mid-span, so each half is a cantilever of length L/2. A simply supported
beam of the same length, loaded at mid-span, has the same two halves: symmetry supplies the
clamp (zero slope and zero slip between the two helix layers at the centre), and the free
ends correspond to the blade's inflection point. With central load P and mid-span sag s, the
stiffness of two blades of length L is

    k = (P/2) / s          (Euler-Bernoulli: 24 B / L^3)

The blade is 128 bp long so that the 59-bp test span has 34 bp of continuing lattice on
each side; a blade with free ends right at the supports splays apart.

## What the models leave out

- No scaffold routing. Each helix has its own scaffold strand, so helix ends are not tied
  by scaffold loops. Real bundle ends are better tied than these.
- oxDNA has no water, no explicit ions and no azobenzene. Contact numbers are a guide only.
- The flexure is not simulated in the contact models. Sideways guidance is imposed with
  harmonic traps, either on the far halves only (`far`) or on every nucleotide more than
  six base pairs from the contact (`ideal`).
- Runs are short (10^5 steps for the bundles). Stiffness from variance is good to about 30%.

## What the full-device check still has to measure

Once the device is drawn in caDNAno, convert it with tacoxDNA, relax it (oxView's
rigid-body relaxation first, for a design with helices in two directions), and run it on
a GPU for about 10^8 steps. Three numbers decide whether staples are ordered:

1. Preload. Simulate the device with and without the stop. The stage shift between the
   two, times the blade stiffness, is the preload. Target 50-80 pN in total and under
   20 pN on each helix end of the stop.
2. Contact jitter. The spread of the distance between the stage face and the stop,
   measured over reference slabs a few base pairs back from each face. Target under 0.15 nm.
3. Step. Simulate stops that are 0, 1 and 2 bp taller. The stage should move by
   0.34 nm per base pair to within 0.1 nm. A smaller step means the stop is sinking.

Run each for long enough to see a collapse if there is one: the single-helix contact here
looked stable for 4 x 10^5 steps at 70 pN and then failed.
