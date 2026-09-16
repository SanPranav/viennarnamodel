# Exon-Skipping Detection & Hairpin Catalytic Assembly (CHA) Modeling

This repository provides a purely ViennaRNA-based Python pipeline to analyze, model, and screen Hairpin Catalytic Assembly (CHA) systems for target detection across exon-skipping splice junctions.

---

## Features

* **Hairpin Folding & Ensemble Analysis:** Computes Minimum Free Energy (MFE) structures, ensemble energies, and ensemble defects/variances for individual hairpins ($H_1$ and $H_2$).
* **Leakage Prediction:** Evaluates 2-strand equilibrium co-folding ($H_1 + H_2$) to assess unintentional background assembly in the absence of a target.
* **Target-Triggered Assembly:** Models 3-strand co-folding ($Target + H_1 + H_2$) to confirm thermodynamically favored, target-catalyzed complex formation.
* **Specificity Screening:** Conducts 3-strand controls with scrambled/mismatched target sequences to verify resistance against off-target triggers.
* **Candidate Sequence Screening:** Evaluates self-folding MFE vs. target duplex binding energy across multiple mRNA candidate targets.

---

## Prerequisites & Installation

### Requirements
* Python 3.8+
* `ViennaRNA` Python library

### Setup Environment

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd <repository-folder>

```

2. **Set up a virtual environment:**
```bash
python3 -m venv .venv
source .venv/bin/activate

```


3. **Install ViennaRNA Python bindings:**
* **via Conda (recommended):**
```bash
conda install -c bioconda viennarna

```


* **via Pip (if available for your platform):**
```bash
pip install viennarna

```





---

## Usage

1. **Configure Parameters:** Update your sequences, temperature, and candidate list in `config.json`.
2. **Execute Analysis:**
```bash
python3 run_analysis.py
```
3. **Interpret Results**
```bash
pyhton3 interpret_results.py
```

---

## File Structure

* `config.json` — System sequence definitions, experimental conditions ($37^\circ\text{C}$ temperature), and candidate target pools.
* `run_analysis.py` — Core evaluation script implementing ViennaRNA multi-strand cofolding (`RNA.fold_compound`, `RNA.duplexfold`, and partition functions).
* `README.md` — Project documentation.
* `interpret_results.py` — Interpreter which compares results through individual hairpins, leak predictions, target-triggered assembly, specificity check.