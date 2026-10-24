# Results of the component checks (4-5 October 2026)

All runs: oxDNA2, 1 M salt, average-sequence parameters, 25 C, dt = 0.005.
Errors are block standard errors. Stiffness from variance is good to about 30%.

## 1. Guide blade: three-point bending

2x8 square-lattice blade, 128 bp long; supports at base-pair planes 34 and 93
(span 59 bp = 19.6 nm), load on planes 63-64; 120,000 steps per run, first 10,000 dropped.

| Central load | Mid-span sag (nm) |
| --- | --- |
| +70 pN | +0.714 +- 0.038 |
| -70 pN | -0.824 +- 0.081 |

- Response: 0.769 +- 0.045 nm per 70 pN, so two blades of this span give
  **0.0455 +- 0.003 N/m** (0.043 N/m scaled to exactly 20 nm).
- Beam-theory estimate in the brief: 0.04-0.07 N/m. The simulation sits at its low end.
- From the thermal spread of the sag alone: 0.069 and 0.039 N/m in the two runs.
- A 70 pN preload needs 1.54 nm of interference (4.5 bp) and stores 13 kT.
- Helix spacing in the span settled within 10,000 steps: 2.7 nm on average, up to 4 nm
  where helices bow between crossovers.
- A 64-bp blade with free ends at the supports splayed apart (up to 8-11 nm between
  neighbouring helix ends) and was not used.

## 2. Stop contact: 3x3 post against 3x3 stage, 70 pN

Ideal sideways guidance (every nucleotide more than 6 bp from the contact held in x and y);
100,000 steps per run, first 10,000 dropped. d is the distance between reference slabs
of the stage and the post.

| Stop | d (nm) | Step above flush (nm) | Jitter in d (nm) | Stiffness (N/m) |
| --- | --- | --- | --- | --- |
| Flush face, nine helix ends | 6.453 +- 0.004 | 0 | 0.040 | 2.6 |
| One helix 2 bp proud | 6.737 +- 0.024 | 0.284 +- 0.024 | 0.113 | 0.32 |
| Four-helix pad 2 bp proud | 7.124 +- 0.014 | 0.670 +- 0.014 | 0.078 | 0.67 |

- Designed step for 2 bp in oxDNA: 0.664 nm. The single helix delivers 43% of it because
  it sinks into its bundle; the four-helix pad delivers all of it.
- The stack itself (gap between the two end base pairs) is stiff in every case: 2-2.6 N/m,
  gap 0.36-0.39 nm. The softness is in how the loaded helix is held by its neighbours.
- With guidance on the far halves only (soft traps, 24 bp from each far end), the
  single-helix contact slid sideways within 40,000 steps and the two faces interleaved.
- No end base pair opened in any of the three guided runs, but see the next section:
  100,000 steps is too short to rule it out for the single helix.

## 3. Smallest contact: two coaxial 16-bp duplexes pressed end to end

Sideways guidance up to 3 bp from the contact; 2,000,000 steps per run.

| Force | Outcome |
| --- | --- |
| 0 pN | stack opened at 470,000 steps; the duplexes drifted apart |
| 18 pN | held |
| 35 pN | held |
| 50 pN | run 1 held; run 2 lost one base pair at about 700,000 steps (sank 0.38 nm) |
| 70 pN | run 1 collapsed at 390,000 steps (sank 1.7 nm); run 2 held |
| 140 pN | collapsed at once |
| 210 pN | crushed at once (all eight watched base pairs open) |

- Reading: the contact needs some preload to stay closed, and fails by an activated
  process whose rate rises with force. The safe range for one helix end in this model is
  roughly 18-35 pN.
- This matches the estimate in `flexure_estimates_check.py` (B6): one base pair's pairing
  free energy spread over its 0.34 nm rise resists about 35 pN at 25 C.
- With 8 bp free on each side instead of 3, the ends slid past each other within
  200,000 steps at 70 pN (`duplex/slip8`).
- oxDNA refuses temperatures below 0 C unless forced (`T_force_value = true`), so the
  planned -40 C run was not done. Cold behaviour has to come from experiment.

## What follows for the design

1. Blade stiffness and preload geometry in the brief stand.
2. The stop should be a pad of about four helix ends, the same ones in both positions,
   with under 20 pN on each. A single-helix piston under 70 pN gives less than half the
   designed step and can collapse.
3. Sideways guidance of the stage at the contact matters as much as the preload. It was
   imposed here; in the device it has to come from the blades, and that is the first
   question for the full-device simulation.
