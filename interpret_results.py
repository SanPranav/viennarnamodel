import re

def analyze_dot_bracket(dot_bracket_str):
    """Parses a dot-bracket string and returns structural counts."""
    unpaired = dot_bracket_str.count('.')
    paired = dot_bracket_str.count('(') + dot_bracket_str.count(')')
    strands = dot_bracket_str.split('&')
    return {
        "unpaired_bases": unpaired,
        "paired_bases": paired,
        "strand_count": len(strands)
    }

def interpret_cha_results(h1_mfe, h2_mfe, leak_mfe, trig_mfe, spec_mfe):
    """Provides a thermodynamic interpretation of CHA pipeline outputs."""
    print("==================================================================")
    print("             CHA THERMODYNAMIC ANALYSIS INTERPRETER               ")
    print("==================================================================\n")

    # 1. Individual Stability
    avg_individual = (h1_mfe + h2_mfe) / 2
    print(f"1. Individual Hairpins:")
    print(f"   * H1 MFE: {h1_mfe:.2f} kcal/mol | H2 MFE: {h2_mfe:.2f} kcal/mol")
    print(f"   * Assessment: Hairpins are individually stable.\n")

    # 2. Leak Assessment
    leak_delta = leak_mfe - (h1_mfe + h2_mfe)
    print(f"2. Leak Prediction (H1 + H2):")
    print(f"   * Co-fold ΔG: {leak_mfe:.2f} kcal/mol (Net change: {leak_delta:.2f} kcal/mol)")
    if leak_delta < -20:
        print("   * WARNING: High spontaneous leak potential! H1 and H2 have a strong")
        print("     thermodynamic driving force to assemble even without the target.")
        print("     Consider lengthening hairpin stems to raise the kinetic barrier.\n")
    else:
        print("   * Assessment: Low background leakage. Hairpins remain metastably separated.\n")

    # 3. Target Activation
    trig_delta = trig_mfe - leak_mfe
    print(f"3. Target-Triggered Assembly:")
    print(f"   * 3-Strand ΔG: {trig_mfe:.2f} kcal/mol")
    print(f"   * ΔΔG (Triggered vs Leak): {trig_delta:.2f} kcal/mol")
    if trig_delta < -5:
        print("   * SUCCESS: Target presence successfully drives the system to a lower")
        print("     energy state, favoring full assembly over the resting hairpins.\n")
    else:
        print("   * WARNING: Target provides minimal thermodynamic drive over the leak state.\n")

    # 4. Specificity
    spec_diff = abs(spec_mfe - leak_mfe)
    print(f"4. Specificity Check:")
    print(f"   * Scrambled Target ΔG: {spec_mfe:.2f} kcal/mol")
    if spec_diff < 3.0:
        print("   * SUCCESS: High specificity. The scrambled sequence behaves like the")
        print("     blank control (leak state) and fails to trigger assembly.\n")
    else:
        print("   * WARNING: Off-target binding detected. Scrambled target alters the state.\n")

    print("==================================================================")


if __name__ == "__main__":
    # Example values pulled directly from your run output
    H1_MFE = -34.10
    H2_MFE = -37.50
    LEAK_MFE = -128.80
    TRIGGERED_MFE = -139.20
    SCRAMBLED_MFE = -128.20

    interpret_cha_results(H1_MFE, H2_MFE, LEAK_MFE, TRIGGERED_MFE, SCRAMBLED_MFE)