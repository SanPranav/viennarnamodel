import json
import os
import sys
import RNA

def load_config(config_file="config.json"):
    """Reads configuration directly from disk."""
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Missing configuration file: {config_file}")
    with open(config_file, "r") as f:
        return json.load(f)

def interpret_cha_results(h1_mfe, h2_mfe, leak_mfe, trig_mfe, spec_mfe):
    """Provides a thermodynamic interpretation using only live-calculated inputs."""
    print("==================================================================")
    print("             CHA THERMODYNAMIC ANALYSIS INTERPRETER               ")
    print("==================================================================\n")

    # 1. Individual Stability
    print("1. Individual Hairpins:")
    print(f"   * H1 MFE: {h1_mfe:.2f} kcal/mol | H2 MFE: {h2_mfe:.2f} kcal/mol")
    print("   * Assessment: Hairpins are individually folded.\n")

    # 2. Leak Assessment
    leak_delta = leak_mfe - (h1_mfe + h2_mfe)
    print("2. Leak Prediction (H1 + H2):")
    print(f"   * Co-fold ΔG: {leak_mfe:.2f} kcal/mol (Net change: {leak_delta:.2f} kcal/mol)")
    if leak_delta < -20:
        print("   * WARNING: High spontaneous leak potential! Strong driving force to assemble without target.\n")
    else:
        print("   * Assessment: Low background leakage. Hairpins remain metastably separated.\n")

    # 3. Target Activation
    trig_delta = trig_mfe - leak_mfe
    print("3. Target-Triggered Assembly:")
    print(f"   * 3-Strand ΔG: {trig_mfe:.2f} kcal/mol")
    print(f"   * ΔΔG (Triggered vs Leak): {trig_delta:.2f} kcal/mol")
    if trig_delta < -5:
        print("   * SUCCESS: Target presence successfully drives the system to a lower energy state.\n")
    else:
        print("   * WARNING: Target provides minimal thermodynamic drive over the leak state.\n")

    # 4. Specificity
    spec_diff = abs(spec_mfe - leak_mfe)
    print("4. Specificity Check:")
    print(f"   * Scrambled Target ΔG: {spec_mfe:.2f} kcal/mol")
    if spec_diff < 3.0:
        print("   * SUCCESS: High specificity. Scrambled sequence fails to trigger assembly.\n")
    else:
        print("   * WARNING: Off-target binding detected. Scrambled target alters the state.\n")

    print("==================================================================")

def run_dynamic_pipeline():
    cfg = load_config("config.json")
    exp = cfg["experimental_conditions"]
    hp = cfg["hairpin_system"]

    # Model settings
    md = RNA.md()
    md.temperature = float(exp["celsius"])

    # Extract sequences live
    seq_h1 = hp["H1_sequence"]
    seq_h2 = hp["H2_sequence"]
    seq_target = hp["target_sequence"]
    seq_scrambled = hp["scrambled_target_sequence"]

    # Calculate MFE live with ViennaRNA
    _, h1_mfe = RNA.fold_compound(seq_h1, md).mfe()
    _, h2_mfe = RNA.fold_compound(seq_h2, md).mfe()
    _, leak_mfe = RNA.fold_compound(f"{seq_h1}&{seq_h2}", md).mfe()
    _, trig_mfe = RNA.fold_compound(f"{seq_target}&{seq_h1}&{seq_h2}", md).mfe()
    _, spec_mfe = RNA.fold_compound(f"{seq_scrambled}&{seq_h1}&{seq_h2}", md).mfe()

    # Pass live-calculated metrics straight to interpreter
    interpret_cha_results(
        h1_mfe=h1_mfe,
        h2_mfe=h2_mfe,
        leak_mfe=leak_mfe,
        trig_mfe=trig_mfe,
        spec_mfe=spec_mfe
    )

if __name__ == "__main__":
    run_dynamic_pipeline()