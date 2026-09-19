import json
import sys
import RNA

# 1. LOAD CONFIGURATION
try:
    with open("config.json", "r") as f:
        config = json.load(f)
except FileNotFoundError:
    print("Error: config.json not found!")
    sys.exit(1)

env_cfg = config["experimental_conditions"]
hp_cfg = config["hairpin_system"]

# Set experimental temperature
RNA.cvar.temperature = env_cfg["celsius"]

# Load clean RNA sequences (convert T to U if present)
seq_h1 = hp_cfg["H1_sequence"].replace("T", "U")
seq_h2 = hp_cfg["H2_sequence"].replace("T", "U")
seq_target = hp_cfg["target_sequence"].replace("T", "U")
seq_luc = hp_cfg["luc_target_sequence"].replace("T", "U")
seq_scrambled = hp_cfg["scrambled_target_sequence"].replace("T", "U")

print("==================================================================")
print("             CHA CONTROLS VALIDATION (ViennaRNA v3)               ")
print("==================================================================\n")

# --- CONTROL 1: Native Target Assembly ---
complex_target = f"{seq_target}&{seq_h1}&{seq_h2}"
fc_target = RNA.fold_compound(complex_target)
struct_target, mfe_target = fc_target.mfe()

# --- CONTROL 2: Luc Target Assembly ---
complex_luc = f"{seq_luc}&{seq_h1}&{seq_h2}"
fc_luc = RNA.fold_compound(complex_luc)
struct_luc, mfe_luc = fc_luc.mfe()

# --- CONTROL 3: Scrambled Target (Off-Target) ---
complex_scram = f"{seq_scrambled}&{seq_h1}&{seq_h2}"
fc_scram = RNA.fold_compound(complex_scram)
struct_scram, mfe_scram = fc_scram.mfe()

# --- CONTROL 4: Un-triggered Leak Baseline ---
complex_leak = f"{seq_h1}&{seq_h2}"
fc_leak = RNA.fold_compound(complex_leak)
struct_leak, mfe_leak = fc_leak.mfe()

# Calculate Net Driving Energies
ddG_target = mfe_target - mfe_leak
ddG_luc = mfe_luc - mfe_leak
ddG_scram = mfe_scram - mfe_leak

print(f"{'Control Condition':<25} | {'3-Strand MFE (kcal/mol)':<23} | {'ΔΔG vs Leak (kcal/mol)':<22}")
print("-" * 75)
print(f"{'Native D2 Target':<25} | {mfe_target:<23.2f} | {ddG_target:<22.2f}")
print(f"{'Luciferase (Luc) Target':<25} | {mfe_luc:<23.2f} | {ddG_luc:<22.2f}")
print(f"{'Scrambled Target (Control)':<25} | {mfe_scram:<23.2f} | {ddG_scram:<22.2f}")
print(f"{'H1 + H2 (Leak Baseline)':<25} | {mfe_leak:<23.2f} | {'0.00 (Baseline)':<22}")
print("-" * 75)

print("\n=== VERIFICATION SUMMARY ===")
if ddG_target < -10.0 and ddG_luc < -10.0:
    print("✓ SUCCESS: Both Native D2 and Luc targets exhibit strong thermodynamic driving force.")
else:
    print("✗ WARNING: Target driving force is weak.")

if mfe_scram > mfe_leak:
    print("✓ SUCCESS: Scrambled control is disfavored relative to baseline (High Specificity).")
else:
    print("✗ WARNING: Off-target sequence triggered false-positive assembly.")

print("\n==================================================================")