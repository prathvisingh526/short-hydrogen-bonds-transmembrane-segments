# Self-contact hydrogen bonds in protein structures

This repository contains the analysis workflow used to identify and characterize
self-contacting hydrogen bonds involving cysteine (Cys), serine (Ser), and
threonine (Thr) residues in high-resolution protein crystal structures, with a
focus on transmembrane segments.

The workflow combines structural data collection, hydrogen-bond detection,
secondary-structure and solvent-accessibility annotation, post-processing, and
preparation of selected systems for quantum-chemical calculations.

## Repository contents

| File | Purpose |
| --- | --- |
| `Short hydrogen bonds in TM segments.ipynb` | Main, documented end-to-end analysis notebook |
| `H-bonding parameters.json` | Donor, hydrogen, acceptor, antecedent, and residue definitions used by the workflow |
| `post_processing.py` | Filtering, classification, summary, and structural post-processing functions |
| `range_to_table.py` | Excel range-to-table formatting helper |
| `residue_counter.py` | Counts relevant residues in transmembrane segments |
| `ssec_strc_phi_psi_chi1_sasa_value_returner.py` | Retrieves secondary structure, backbone/side-chain angles, and solvent-accessibility values |
| `chimerax_phi_psi_chi1_files_generator.py` | UCSF ChimeraX script for exporting phi, psi, and chi1 attributes |
| `chimera_hadding_script_for_short_hydrogen_bonds.py` | UCSF Chimera script used to add hydrogens to selected structures |

## Python setup

Python 3.12 was recorded in the notebook metadata. Create an isolated
environment and install the Python dependencies:

```bash
python -m venv .venv
```

On Windows:

```powershell
.venv\Scripts\activate
python -m pip install -r requirements.txt
jupyter lab
```

On Linux or macOS:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
jupyter lab
```

Open `Short hydrogen bonds in TM segments.ipynb` and run the sections required
for your analysis.

## External scientific software

Some stages are not ordinary Python processes and require separate
installations:

- UCSF Chimera and/or UCSF ChimeraX
- MolProbity `Reduce`
- STRIDE
- NACCESS
- Gaussian, for the quantum-chemical stages
- WSL with Ubuntu, where used by the original Windows workflow

These programs are not installed by `requirements.txt`; follow their official
installation and licensing instructions.

## Important portability note

The notebook and `ssec_strc_phi_psi_chi1_sasa_value_returner.py` retain the
original absolute Windows paths (for example, paths on `C:`, `D:`, and `G:`)
and executable locations. Before running the workflow on another computer,
search for these path assignments and replace them with locations appropriate
to your datasets and software installation.

The notebook represents a research workflow rather than a single one-click
pipeline. Several sections are explicitly manual, depend on outputs from prior
sections, or call external databases and desktop scientific applications.
Run it sequentially while checking the explanatory Markdown cells.

## Data sources and network access

The notebook contains stages that retrieve or process information from
structural biology resources including RCSB PDB, PDBTM, OPM, and UniProt.
Availability and response formats can change, so verify the relevant service
documentation before rerunning older retrieval cells.

## Reproducibility

- Keep the notebook and helper modules in the same directory so that local
  imports continue to work.
- Preserve `H-bonding parameters.json` alongside the notebook or update the
  corresponding path in the relevant cell.
- Large downloaded structures, generated spreadsheets, logs, and Gaussian
  outputs are intentionally excluded from version control by `.gitignore`.
- Record software versions and input-data download dates when producing a new
  analysis run.

## Validation performed for this repository

- The notebook and parameter file were checked as valid JSON.
- All included Python files were checked for Python syntax errors.
- Local module imports referenced by the notebook are included in the
  repository.

Full execution was not attempted because the original input datasets,
proprietary/external scientific programs, and original absolute paths are not
included.

## Citation and license

If this workflow supports a publication, add the corresponding paper and/or
dataset citation here before making the repository public. No open-source
license has been selected; add a `LICENSE` file only after choosing terms that
match the authors' and institution's requirements.
