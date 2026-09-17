import RNA
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. SEQUENCES & EXPERIMENTAL CONDITIONS (FROM CONFIG)
# ------------------------------------------------------------------------------
temp_celsius = 37.0
H1_0 = 1e-07          # Molar
H2_0 = 5e-08          # Molar
Target_0 = 2.5e-08      # Molar

# Hairpins and Candidate 28 Target
H1_seq = "GCUACCAGGUAGCCCCAGCCGCCACGUCUUGCGGCUGGGGCUACCUGCGAGACAACAGAGAGACGGUCGGGUCCACAGUUUCCGGAAACUGUGU"
H2_seq = "CGAGUAGAGUGUGGGCUCUCUGUUGUCUCGAGAUGUCGCAAGACGUGGCGGCUGGGGCUACCUGCCAUGUCUUGGACAUCUACCAACAGUAUCGCUC"

# Test Targets
target_seq = "GCUACCAGGUAGCCCCAGC"                           # Candidate 28 Target
intron_junction_seq = "ACAGGUAAGUAUCAAGGUUACAAG"                  # Exon1-Intron Pre-mRNA Junction (bp 2340)[cite: 1]
scrambled_seq = "UUUUUUUUUUUUUUUUUUU"                         # Scrambled Control

# Set ViennaRNA global temperature
RNA.cvar.temperature = temp_celsius

# ------------------------------------------------------------------------------
# 2. THERMODYNAMIC EVALUATION (ViennaRNA)
# ------------------------------------------------------------------------------
def evaluate_binding(target, hairpin):
    duplex = f"{target}&{hairpin}"
    fc = RNA.fold_compound(duplex)
    _, mfe = fc.mfe()
    return mfe

dg_target = evaluate_binding(target_seq, H1_seq)
dg_intron = evaluate_binding(intron_junction_seq, H1_seq)
dg_scrambled = evaluate_binding(scrambled_seq, H1_seq)

print(f"=== THERMODYNAMIC BINDING ENERGIES (ΔG at {temp_celsius}°C) ===")
print(f"Target Binding ΔG           : {dg_target:.2f} kcal/mol")
print(f"Intron Junction Binding ΔG  : {dg_intron:.2f} kcal/mol")
print(f"Scrambled Control Binding ΔG: {dg_scrambled:.2f} kcal/mol\n")

# ------------------------------------------------------------------------------
# 3. KINETIC SIMULATION COMPARISON
# ------------------------------------------------------------------------------
# Kinetic rates based on duplex stability relative to target
k1_target = 1e5    # Target bimolecular rate (M^-1 s^-1)
k1_intron = 1e2    # Off-target pre-mRNA binding rate
k1_scram = 10      # Scrambled non-specific rate
k1_r = 0.1         # Unbinding rate (s^-1)
k2 = 1e5           # H2 displacement rate (M^-1 s^-1)
k_leak = 10        # Spontaneous leakage rate (M^-1 s^-1)

def simulate_cha(k1_rate, t_span=(0, 7200)):
    def cha_system(t, y):
        T, H1, H2, H1T, C3, Leak = [max(0.0, val) for val in y]
        dT = -k1_rate * T * H1 + k1_r * H1T + k2 * H1T * H2
        dH1 = -k1_rate * T * H1 + k1_r * H1T - k_leak * H1 * H2
        dH2 = -k2 * H1T * H2 - k_leak * H1 * H2
        dH1T = k1_rate * T * H1 - k1_r * H1T - k2 * H1T * H2
        dC3 = k2 * H1T * H2
        dLeak = k_leak * H1 * H2
        return [dT, dH1, dH2, dH1T, dC3, dLeak]

    y0 = [Target_0, H1_0, H2_0, 0, 0, 0]
    t_eval = np.linspace(t_span[0], t_span[1], 500)
    sol = solve_ivp(cha_system, t_span, y0, t_eval=t_eval, method='Radau')
    
    signal = sol.y[4] * 1e9  # nM
    leak = sol.y[5] * 1e9    # nM
    return sol.t / 60, signal, leak

t_min, sig_target, leak_bg = simulate_cha(k1_target)
_, sig_intron, _ = simulate_cha(k1_intron)
_, sig_scram, _ = simulate_cha(k1_scram)

# ------------------------------------------------------------------------------
# 4. PLOTTING & OUTPUT
# ------------------------------------------------------------------------------
plt.figure(figsize=(9, 5))
plt.plot(t_min, sig_target, label="Candidate 28 Target", color="green", linewidth=2)
plt.plot(t_min, sig_intron, label="Pre-mRNA Intron Junction", color="orange", linestyle="-.", linewidth=2)
plt.plot(t_min, sig_scram, label="Scrambled (U19)", color="blue", linestyle=":", linewidth=2)
plt.plot(t_min, leak_bg, label="Spontaneous Leak Baseline", color="red", linestyle="--", linewidth=1.5)

plt.xlabel("Time (minutes)")
plt.ylabel("Active Complex Concentration (nM)")
plt.title("CHA Kinetics: Candidate 28 vs. Junction Controls (37°C)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("junction_controls_kinetics.png")

print("=== KINETIC SIMULATION RESULTS (2 Hours) ===")
print(f"Target Signal     : {sig_target[-1]:.2f} nM")
print(f"Intron Signal     : {sig_intron[-1]:.2f} nM")
print(f"Scrambled Signal  : {sig_scram[-1]:.2f} nM")
print(f"Leak Baseline     : {leak_bg[-1]:.2f} nM")
print("\nPlot saved as 'junction_controls_kinetics.png'.")