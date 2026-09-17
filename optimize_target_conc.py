import json
import numpy as np
import RNA

# Load system configuration
with open("config.json", "r") as f:
    config = json.load(f)

hp_cfg = config["hairpin_system"]
env_cfg = config["experimental_conditions"]

RNA.cvar.temperature = env_cfg["celsius"]

seq_h1 = hp_cfg["H1_sequence"].replace("T", "U")
seq_h2 = hp_cfg["H2_sequence"].replace("T", "U")
seq_target = hp_cfg["target_sequence"].replace("T", "U")

fixed_h1_m = env_cfg["concentration_H1_M"]
fixed_h2_m = env_cfg["concentration_H2_M"]

# Range of target concentrations to test (10 pM to 1 uM)
target_concentrations = np.logspace(-11, -6, num=10)

def compute_complex_deltag(seq_t, seq_1, seq_2):
    complex_seq = f"{seq_t}&{seq_1}&{seq_2}"
    fc = RNA.fold_compound(complex_seq)
    _, mfe = fc.mfe()
    return mfe

mfe_complex = compute_complex_deltag(seq_target, seq_h1, seq_h2)

print("==================================================================")
print("     TARGET CONCENTRATION OPTIMIZATION SWEEP                     ")
print("==================================================================")
print(f"Fixed H1 Concentration     : {fixed_h1_m:.1e} M")
print(f"Fixed H2 Concentration     : {fixed_h2_m:.1e} M")
print(f"3-Strand Assembly MFE ΔG   : {mfe_complex:.2f} kcal/mol\n")

print(f"{'Target Conc (M)':<18} | {'Target : H1 Ratio':<18} | {'Regime Status':<20}")
print("-" * 62)

for target_m in target_concentrations:
    ratio = target_m / fixed_h1_m
    
    if ratio < 0.1:
        status = "Sub-stoichiometric (Low)"
    elif 0.1 <= ratio <= 1.0:
        status = "Optimal Catalytic Range"
    elif 1.0 < ratio <= 2.0:
        status = "Equimolar Saturation"
    else:
        status = "Target Excess (Sub-optimal)"
        
    print(f"{target_m:<18.2e} | {ratio:<18.2f} | {status:<20}")

print("\n==================================================================")