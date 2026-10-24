"""
flexure_estimates_check.py
Independent rebuild of the estimates in "Single-Axis DNA Nanopositioner: Experiment Brief".

Written without sight of flexure_estimates.py, from the numbers in the brief, so that the
two can be compared. Part A reproduces the brief. Part B adds three effects that the
beam-theory / hard-wall treatment leaves out.

Units: SI internally; printed in pN, nm, N/m.
"""
import math

kB = 1.380649e-23
T_HOT, T_COLD = 298.15, 233.15          # 25 C, -40 C
nm, pN = 1e-9, 1e-12

# ---------------------------------------------------------------- geometry and materials
R_HELIX = 1.0 * nm                      # mechanical radius of one duplex
BP_RISE = 0.34 * nm
N_THICK, N_WIDE = 2, 8                  # blade cross-section, helices
L_BLADE = 20 * nm
N_BLADES = 2
D_SQ = 2.5 * nm                         # inter-helix spacing, square lattice (2.25-2.6 reported)
D_HC = 2.25 * nm                        # inter-helix spacing, honeycomb (for the 6HB reference)
LP_6HB = (1.88e-6, 3.3e-6)              # six-helix-bundle persistence length: measured low, rigid-coupling high
LP_DUPLEX = 50 * nm
S_DUPLEX = 1100 * pN                    # stretch modulus of one duplex
F_OVERSTRETCH = 65 * pN
PRELOAD = 70 * pN
TARGET = 0.1 * nm
R0 = 5.5 * nm

A_H = math.pi * R_HELIX**2
I_H = math.pi * R_HELIX**4 / 4


def kT(T):
    return kB * T


def I_6hb(d=D_HC):
    # six helices on a hexagon of circumradius d: sum(y^2) = 3 d^2
    return 6 * I_H + A_H * 3 * d**2


def I_blade(d=D_SQ, coupled=True):
    # bending about the wide axis; the two layers sit at +-d/2 from the neutral plane
    n = N_THICK * N_WIDE
    return n * I_H + (n * A_H * (d / 2) ** 2 if coupled else 0.0)


def B_blade(lp_6hb, T=T_HOT, coupled=True):
    # "beam theory calibrated to six-helix-bundle stiffness": scale B by the ratio of second moments
    B6 = kT(T_HOT) * lp_6hb             # calibrated at room temperature
    return B6 * I_blade(coupled=coupled) / I_6hb()


def k_parallelogram(B, L=L_BLADE, n=N_BLADES, shear_rigidity=None):
    """Two fixed-guided blades. Timoshenko correction if a shear rigidity (N) is given."""
    k = 12 * B / L**3
    if shear_rigidity is not None:
        k /= 1 + 12 * B / (shear_rigidity * L**2)
    return n * k


def line(label, value, unit=""):
    print(f"  {label:<62s} {value:>10s} {unit}")


print("=" * 86)
print("PART A  -  reproducing the brief")
print("=" * 86)

print("\nA1. Cooling and thermal jitter")
ratio = math.sqrt(T_COLD / T_HOT)
line("jitter(-40 C) / jitter(25 C) = sqrt(T2/T1)", f"{ratio:.3f}")
line("reduction", f"{100*(1-ratio):.1f}", "%   (brief: about 12%)")

print("\nA2. Preload for mean lift-off < 0.1 nm against a hard stop:  F >= kT/delta")
for T, lab in ((T_HOT, "25 C"), (T_COLD, "-40 C")):
    line(f"F_min at {lab}", f"{kT(T)/TARGET/pN:.1f}", "pN   (brief: 41 / 32)")

print("\nA3. Stiffness for rms jitter < 0.1 nm:  k >= kT/sigma^2")
for T, lab in ((T_HOT, "25 C"), (T_COLD, "-40 C")):
    line(f"k_min at {lab}", f"{kT(T)/TARGET**2:.3f}", "N/m  (brief: 0.41 at 25 C)")

print("\nA4. Blade stiffness along the motion, Euler-Bernoulli, layers rigidly coupled")
line("I(6HB), honeycomb d = 2.25 nm", f"{I_6hb()/nm**4:.1f}", "nm^4")
line("I(2x8 blade), square lattice d = 2.5 nm", f"{I_blade()/nm**4:.1f}", "nm^4")
line("ratio", f"{I_blade()/I_6hb():.2f}")
ks = []
for lp in LP_6HB:
    B = B_blade(lp)
    k = k_parallelogram(B)
    ks.append(k)
    line(f"Lp(6HB) = {lp*1e6:.2f} um -> B_blade", f"{B:.2e}", "N m^2")
    line("   two-blade stiffness 24 B / L^3", f"{k:.3f}", "N/m  (brief: 0.04-0.07)")
k_lo, k_hi = min(ks), max(ks)

print("\nA5. Interference needed for 70 pN preload")
for k in (k_lo, k_hi):
    d = PRELOAD / k
    line(f"k = {k:.3f} N/m -> deflection", f"{d/nm:.2f}", f"nm = {d/BP_RISE:.1f} bp  (brief: 3-6 bp)")
    line("   strain energy stored, F*delta/2", f"{0.5*PRELOAD*d/kT(T_HOT):.1f}", "kT")

print("\nA6. Peak axial force in a blade helix at 70 pN")
M = (PRELOAD / N_BLADES) * L_BLADE / 2                       # root moment, fixed-guided
f_helix = M * (D_SQ / 2) * A_H / I_blade()
line("root moment per blade", f"{M/(pN*nm):.0f}", "pN nm")
line("force per helix", f"{f_helix/pN:.1f}", "pN   (brief: about 15; overstretch about 65)")

print("\nA7. FRET change for a one-base-pair step at r = R0")
dEdr = 1.5 / R0                                               # |dE/dr| at r = R0
line("dE per 0.34 nm, R0 = 5.5 nm", f"{dEdr*BP_RISE:.3f}", "     (brief: about 0.10)")

print()
print("=" * 86)
print("PART B  -  what the beam-theory / hard-wall treatment leaves out")
print("=" * 86)

print("\nB1. The blade is short: L/h = 20 nm / 5 nm = 4, with about three crossover planes.")
print("    Lower bound if the two layers bend independently (no shear transfer):")
for lp in LP_6HB:
    k = k_parallelogram(B_blade(lp, coupled=False))
    line(f"Lp(6HB) = {lp*1e6:.2f} um -> uncoupled two-blade stiffness", f"{k:.4f}", "N/m")
line("coupled / uncoupled", f"{I_blade()/I_blade(coupled=False):.1f}", "x")
print("    So the plausible range is about 0.006-0.07 N/m, and the interference for 70 pN")
print("    is 1-12 nm (3-35 bp). oxDNA has to fix this number before the post length is set.")

print("\nB2. The stop is not a hard wall. For a repulsion that decays over a length lam,")
print("    the contact stiffness at preload F is F/lam and the jitter is sqrt(kT*lam/F).")
print("    (A hard wall is the limit lam -> kT/F = %.3f nm at 70 pN.)" % (kT(T_HOT) / PRELOAD / nm))
print(f"    {'decay length':<34s}{'k_contact':>10s}{'jitter 25C':>12s}{'jitter -40C':>13s}{'F for 0.1 nm':>14s}")
for lab, lam in (("0.1 nm (steric)", 0.1), ("0.3 nm (hydration, Rau-Parsegian)", 0.3),
                 ("1.0 nm (screened electrostatic)", 1.0)):
    lam_m = lam * nm
    kc = PRELOAD / lam_m
    j1 = math.sqrt(kT(T_HOT) * lam_m / PRELOAD)
    j2 = math.sqrt(kT(T_COLD) * lam_m / PRELOAD)
    Fneed = kT(T_HOT) * lam_m / TARGET**2
    print(f"    {lab:<34s}{kc:>8.2f} N/m{j1/nm:>9.3f} nm{j2/nm:>10.3f} nm{Fneed/pN:>11.0f} pN")
print("    A base-stacked (coaxial) contact is stiff by comparison: one stacking step")
kstack = S_DUPLEX / BP_RISE
line("   k = S / rise", f"{kstack:.1f}", "N/m")
line("   jitter at 25 C", f"{math.sqrt(kT(T_HOT)/kstack)/nm:.3f}", "nm")

print("\nB3. The piston carries the whole preload on one duplex.")
for nbp in (1, 2):
    w = PRELOAD * nbp * BP_RISE
    line(f"work done if the stage drops {nbp} bp ({nbp*0.34:.2f} nm)",
         f"{w/kT(T_HOT):.1f}", f"kT at 25 C = {w/kT(T_COLD):.1f} kT at -40 C = {w/(pN*nm)/6.95:.1f} kcal/mol")
print("    A terminal base pair is worth roughly 1-3 kcal/mol, so the protruding pairs are")
print("    not held by their own pairing; the staple body must hold them, and the tip must")
print("    not be able to fray out of the way. Check in simulation.")
print("    Axial load path: stiffness of n bp of single duplex = S / (n * 0.34 nm)")
for n in (4, 7, 10, 16):
    line(f"   n = {n} bp", f"{S_DUPLEX/(n*BP_RISE):.2f}", "N/m")
line("longest unsupported path for 0.41 N/m", f"{S_DUPLEX/0.412/BP_RISE:.1f}", "bp")
print("    Euler buckling of the unsupported duplex (pinned ends):")
for n in (8, 16, 32):
    Lseg = n * BP_RISE
    Fc = math.pi**2 * kT(T_HOT) * LP_DUPLEX / Lseg**2
    line(f"   {n} bp between crossovers", f"{Fc/pN:.0f}", "pN")

print("\nB4. The structural loop. A member of length L that carries the stop sideways off")
print("    the base bends as a cantilever, k = 3B/L^3, in series with the contact.")
for nlay in (2, 4, 8):
    ys = [(i - (nlay - 1) / 2) * D_SQ for i in range(nlay)]
    I = 8 * (nlay * I_H + A_H * sum(y * y for y in ys))
    B = kT(T_HOT) * LP_6HB[0] * I / I_6hb()
    line(f"   {nlay}x8 helices, L = 25 nm", f"{3*B/(25*nm)**3:.3f}", "N/m")
print("    Even 64 helices do not reach 0.41 N/m, so 0.1 nm can only be claimed for the stage")
print("    relative to the stop, not relative to a distant frame. The FRET pair spanning")
print("    the contact measures the right thing; the brief should say so.")

print("\nB5. Readout in two-thirds glycerol")
n_w, n_g = 1.333, 1.42
R0g = R0 * (n_w / n_g) ** (2 / 3)
line("R0 corrected for refractive index 1.42", f"{R0g/nm:.2f}", "nm")
line("dE per 0.34 nm at r = R0", f"{1.5/R0g*BP_RISE:.3f}")
for N in (200, 1000, 5000):
    sE = 0.5 / math.sqrt(N)
    line(f"shot-noise sigma_E at E = 0.5, {N} photons",
         f"{sE:.3f}", f"= {sE/(1.5/R0g)/nm:.3f} nm")

print("\nB6. What one base pair can hold. Pairing free energy per step, spread over its rise,")
print("    is the largest axial load a helix end can carry before its end pairs are squeezed out.")
print("    Nearest-neighbour averages (SantaLucia 1998, 1 M NaCl; values typed from memory, check before citing):")
print("    dH = -8.4 kcal/mol, dS = -22.4 cal/mol/K;")
print("    extrapolated to -40 C with constant dH and dS, which is far outside the measured range.")
KCAL = 6.948                                                  # 1 kcal/mol in pN nm
for lab, dH, dS in (("average step", -8.36, -22.37), ("weakest (TA/AT)", -7.2, -21.3), ("strongest (GC/CG)", -9.8, -24.4)):
    out = []
    for T in (T_HOT, T_COLD):
        dG = dH - T * dS / 1000.0
        out.append("%5.2f kcal/mol = %4.1f pN" % (dG, -dG * KCAL / (BP_RISE / nm)))
    print(f"    {lab:<20s} 25 C: {out[0]}     -40 C: {out[1]}")
print("    So one helix end holds about 35 pN at 25 C on average; 70 pN needs the load shared by")
print("    several helix ends, or a colder sample. oxDNA agrees in kind: two duplexes pressed end to")
print("    end held at 18 and 35 pN, and at 50 and 70 pN half the runs lost base pairs or collapsed")
print("    (see sim/RESULTS.md).")
