<img src="logo.png" alt="Alt Text" width="200">

[Webserver](https://services.healthtech.dtu.dk/services/SIMAlign-1.0/)

SIMalign is a Python command-line tool for protein structure alignment and hotspot prediction. It integrates both automated homology search via Foldseek and user-specified templates, performs structural superposition, calculates similarity scores, and identifies hotspots suitable for mutagenesis. Results are visualized in PyMOL and exported in multiple formats for downstream analysis.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Features

* Automatic homology search with **Foldseek** or user-specified templates
* Supports `.pdb` and `.cif` file formats
* Gap recognition and RMSD filtering for precise alignment
* Similarity scoring using customizable **BLOSUM** matrices
* Hotspot prediction for single and double mutations
* PyMOL session generation (`.pse`) for interactive visualization
* Outputs: JSON scores, HTML hotspot reports, Clustal alignments

## Installation

### 1. Clone the repository
> **Note:** Do **not** include any spaces in the path where you clone the repo.

```bash
git clone https://github.com/morth-lab/SIMalign.git
cd SIMalign
```

### 2. Conda environment

```bash
conda env create -f environment.yml
conda activate simalign_env
```


## Usage

```bash
python scripts/simalign --QUERY query.pdb [options]
```


### Common Options

| Option (`-short`)               | Description                                                        | Default       |
| ------------------------------- | ------------------------------------------------------------------ | ------------- |
| `--QUERY` `-q`                  | Path to input structure file (`.pdb` or `.cif`). **Required**      | —             |
| `--TEMPLATES` `-t`              | One or more template files (for `user_specified`).                 | `None`        |
| `--TEMPLATES_DIR` `-t-dir`      | Directory of template files (for `user_specified`).                | `None`        |
| `--HOMOLOGY_SEARCH_METHOD` `-H` | `foldseek` or `user_specified`.                                    | `foldseek`    |
| `--MAX_DISTANCE` `-d`           | Distance threshold for gap detection (Å).                          | `7`           |
| `--MAX_RMSD` `-r`               | Maximum RMSD for template filtering (Å).                           | `5`           |
| `--FOLDSEEK_DATABASES` `-fd`    | Foldseek DBs (`afdb50`,`afdb_swissprot`,`afdb_proteome`,`pdb100`). | `afdb50`      |
| `--FOLDSEEK_MODE` `-fm`         | `tmalign` or `3diaa`.                                              | `tmalign`     |
| `--FOLDSEEK_THRESHOLD` `-ft`    | Foldseek score/E-value threshold.                                  | `0.7`         |
| `--NUMB_TEMPLATES` `-nt`        | Number of top templates.                                           | `20`          |
| `--BLOSUM` `-b`                 | BLOSUM matrix (`BLOSUM50`,`BLOSUM62`).                             | `BLOSUM62`    |
| `--RESULT_DIR` `-R`             | Directory for results.                                             | `./<JOB_KEY>` |
| `--TMP_DIR` `-tmp`              | Directory for temporary files.                                     | `./tmp`       |
| `--JOB_KEY` `-j`                | Job name key. Auto-generated if omitted.                           | random        |
| `--only_core`  	                | If set to 1, only hotspots in the core of the protein will be considered.  | `1` |

<!-- | `--SEQUENCE_IDENTITY` `-sident` | Min. identity for BLASTp (0–1).                                    | `0.6`         |
| `--SEQUENCE_COV` `-scov`        | Min. coverage for BLASTp (0–1).                                    | `0.6`         |
| `--E_VALUE` `-e`                | E-value threshold for BLASTp.                                      | `0.001`       |
| `--REDUNDANCY_THRESHOLD` `-rt`  | MSA redundancy threshold (0–1).                                    | `0.9`         | -->

Run `simalign --help` for full details.

## Output

When complete, the output directory contains:

* `alignment.aln` – Multiple sequence alignment (Clustal format)
* `scores.json` – Per-residue similarity scores
* `hotspots_mode_1.html` – Single mutation hotspot report
* `hotspots_mode_2.html` – Double mutation hotspot report
* `SIMalign_<JOB_KEY>.pse` – PyMOL session file

## Citation

Please cite our [publication](https://services.healthtech.dtu.dk/services/SIMAlign-1.0/) if you use SIMalign:

```bibtex
@article{ostergaard2025simplicity,
  title={SIMalign: Structure-based alignment and hotspot prediction for protein engineering},
  author={Ostergaard, M. et al.},
  journal={Journal of Molecular Biology},
  year={2025},
  doi={10.1016/j.jmb.2025.03.012}
}
```

## License

Distributed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

* **PyMOL Script Repository** for `findSurfaceResidues.py`
* Developed with support from ChatGPT (May 2025)
