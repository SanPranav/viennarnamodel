import json
import os
import sys
import RNA

def load_config(config_file="config.json"):
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            return json.load(f)
    raise FileNotFoundError(f"Configuration file {config_file} not found.")

def wallace_tm(seq):
    seq = seq.upper()
    a_t = seq.count('A') + seq.count('T') + seq.count('U')
    g_c = seq.count('G') + seq.count('C')
    return 2 * a_t + 4 * g_c

def evaluate_cha(config):
    exp = config["experimental_conditions"]
    hp = config["hairpin_system"]
    
    print("==========================================================")
    print("           CHA SYSTEM ANALYSIS (JSON CONFIG)              ")
    print("==========================================================")
    print(f"Temperature : {exp['celsius']} °C")
    print(f"H1 Conc     : {exp['concentration_H1_M']} M")
    print(f"H2 Conc     : {exp['concentration_H2_M']} M")
    print(f"Target Conc : {exp['concentration_target_M']} M\n")

    # Set up ViennaRNA model details
    md = RNA.md()
    md.temperature = float(exp["celsius"])

    # Individual Hairpin MFEs
    h1_fc = RNA.fold_compound(hp["H1_sequence"], md)
    _, h1_mfe = h1_fc.mfe()
    
    h2_fc = RNA.fold_compound(hp["H2_sequence"], md)
    _, h2_mfe = h2_fc.mfe()

    # Leak State (H1 + H2)
    leak_seq = f"{hp['H1_sequence']}&{hp['H2_sequence']}"
    leak_fc = RNA.fold_compound(leak_seq, md)
    _, leak_mfe = leak_fc.mfe()

    # Triggered State (Target + H1 + H2)
    trig_seq = f"{hp['target_sequence']}&{hp['H1_sequence']}&{hp['H2_sequence']}"
    trig_fc = RNA.fold_compound(trig_seq, md)
    _, trig_mfe = trig_fc.mfe()

    # Scrambled Control (Scrambled + H1 + H2)
    spec_seq = f"{hp['scrambled_target_sequence']}&{hp['H1_sequence']}&{hp['H2_sequence']}"
    spec_fc = RNA.fold_compound(spec_seq, md)
    _, spec_mfe = spec_fc.mfe()

    print("--- THERMODYNAMIC EVALUATION ---")
    print(f"H1 MFE                  : {h1_mfe:.2f} kcal/mol")
    print(f"H2 MFE                  : {h2_mfe:.2f} kcal/mol")
    print(f"H1+H2 Leak ΔG          : {leak_mfe:.2f} kcal/mol")
    print(f"Triggered 3-Strand ΔG   : {trig_mfe:.2f} kcal/mol")
    print(f"Scrambled 3-Strand ΔG   : {spec_mfe:.2f} kcal/mol")
    print(f"Target ΔΔG (Drive)      : {trig_mfe - leak_mfe:.2f} kcal/mol\n")

def evaluate_candidates(config):
    screen_cfg = config.get("candidate_screening", {})
    
    # Check if we are running a single selected_strand or a list of candidates
    if "selected_strand" in screen_cfg:
        cand = screen_cfg["selected_strand"]
        candidates = [cand]  # Wrap single strand in a list so the loop still works
    elif "candidates" in screen_cfg:
        candidates = screen_cfg["candidates"]
    else:
        print("Error: No 'selected_strand' or 'candidates' key found in config.json")
        return

    target_rna = config["hairpin_system"]["target_sequence"].replace("T", "U")

    for cand in candidates:
        cand_seq = cand["sequence"].replace("T", "U")
        cand_idx = cand.get("index", "N/A")
        
        # Calculate folding energy
        fc = RNA.fold_compound(cand_seq)
        _, mfe = fc.mfe()
        
        print(f"Candidate {cand_idx} ({cand_seq}): Self-MFE = {mfe:.2f} kcal/mol")

if __name__ == "__main__":
    cfg = load_config("config.json")
    evaluate_cha(cfg)
    evaluate_candidates(cfg)