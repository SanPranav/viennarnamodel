import json
import sys
import RNA

# ------------------------------------------------------------------------------
# 1. LOAD CONFIGURATION
# ------------------------------------------------------------------------------
try:
    with open("config.json", "r") as f:
        config = json.load(f)
except FileNotFoundError:
    print("Error: config.json not found!")
    sys.exit(1)

env_cfg = config["experimental_conditions"]
hp_cfg = config["hairpin_system"]
screen_cfg = config["candidate_screening"]

# Configure model details dynamically
md = RNA.md()
md.temperature = float(env_cfg["celsius"])

print("==================================================================")
print("     EXON-SKIPPING DETECTION ANALYSIS (ViennaRNA Full Pipeline)   ")
print("==================================================================\n")

seq_h1 = hp_cfg["H1_sequence"].replace("T", "U")
seq_h2 = hp_cfg["H2_sequence"].replace("T", "U")
seq_target = hp_cfg["target_sequence"].replace("T", "U")
seq_scrambled = hp_cfg["scrambled_target_sequence"].replace("T", "U")

# Helper function for partition function and ensemble energy using dynamic model details
def get_ensemble_info(sequence, model_details):
    fc = RNA.fold_compound(sequence, model_details)
    struct, mfe_energy = fc.mfe()
    _, ensemble_energy = fc.pf()
    ensemble_defect = ensemble_energy - mfe_energy
    return struct, mfe_energy, ensemble_energy, ensemble_defect

# ------------------------------------------------------------------------------
# TASK 1: HAIRPIN FOLDING & ENSEMBLE ANALYSIS
# ------------------------------------------------------------------------------
print("=== TASK 1: Hairpin Folding & Ensemble Analysis ===")

struct_h1, mfe_h1, ens_h1, def_h1 = get_ensemble_info(seq_h1, md)
struct_h2, mfe_h2, ens_h2, def_h2 = get_ensemble_info(seq_h2, md)

print(f"H1 MFE Energy       : {mfe_h1:.2f} kcal/mol")
print(f"H1 Ensemble Energy  : {ens_h1:.2f} kcal/mol (Defect / Variance: {def_h1:.2f})")
print(f"H1 Structure        : {struct_h1}\n")

print(f"H2 MFE Energy       : {mfe_h2:.2f} kcal/mol")
print(f"H2 Ensemble Energy  : {ens_h2:.2f} kcal/mol (Defect / Variance: {def_h2:.2f})")
print(f"H2 Structure        : {struct_h2}\n")

# ------------------------------------------------------------------------------
# TASK 2: LEAK PREDICTION (2-Strand System)
# ------------------------------------------------------------------------------
print("=== TASK 2: Leak Prediction (H1 + H2) ===")
leak_complex = f"{seq_h1}&{seq_h2}"
struct_leak, mfe_leak, ens_leak, _ = get_ensemble_info(leak_complex, md)

print(f"H1:H2 MFE ΔG        : {mfe_leak:.2f} kcal/mol")
print(f"H1:H2 Ensemble ΔG   : {ens_leak:.2f} kcal/mol")
print(f"Structure           : {struct_leak}\n")

# ------------------------------------------------------------------------------
# TASK 3: TARGET-TRIGGERED ASSEMBLY (3-Strand System)
# ------------------------------------------------------------------------------
print("=== TASK 3: Target-Triggered Assembly (Target + H1 + H2) ===")
assembly_complex = f"{seq_target}&{seq_h1}&{seq_h2}"
struct_trig, mfe_trig, ens_trig, _ = get_ensemble_info(assembly_complex, md)

print(f"3-Strand MFE ΔG     : {mfe_trig:.2f} kcal/mol")
print(f"3-Strand Ensemble ΔG: {ens_trig:.2f} kcal/mol")
print(f"Structure           : {struct_trig}\n")

# ------------------------------------------------------------------------------
# TASK 4: SPECIFICITY CHECK (3-Strand System with Scrambled Target)
# ------------------------------------------------------------------------------
print("=== TASK 4: Specificity Check (Scrambled Target + H1 + H2) ===")
spec_complex = f"{seq_scrambled}&{seq_h1}&{seq_h2}"
struct_spec, mfe_spec, ens_spec, _ = get_ensemble_info(spec_complex, md)

print(f"Scrambled MFE ΔG    : {mfe_spec:.2f} kcal/mol")
print(f"Scrambled Ensemble ΔG: {ens_spec:.2f} kcal/mol")
print(f"Structure           : {struct_spec}\n")

# ------------------------------------------------------------------------------
# PART B: VIENNARNA CANDIDATE SCREENING
# ------------------------------------------------------------------------------
print("=== PART B: Candidate Screening ===")

selected = screen_cfg.get("selected_strand", {})
strand_idx = selected.get("index", "N/A")
cand_seq = selected.get("sequence", seq_target).replace("T", "U")

fc_cand = RNA.fold_compound(cand_seq, md)
_, cand_mfe = fc_cand.mfe()

# Calculate duplex energy using fold_compound to preserve temperature settings
duplex_seq = f"{cand_seq}&{seq_h1}"
fc_duplex = RNA.fold_compound(duplex_seq, md)
_, duplex_mfe = fc_duplex.mfe()

print(f"Selected Strand Index : {strand_idx}")
print(f"Strand Sequence       : {cand_seq}")
print(f"Strand Self-MFE       : {cand_mfe:.2f} kcal/mol")
print(f"Duplex ΔG with H1     : {duplex_mfe:.2f} kcal/mol")

print("\n==================================================================")
print("                    ANALYSIS COMPLETE                             ")
print("==================================================================")