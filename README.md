# SIMalign

SIMalign is a Python tool designed for protein structure alignment and analysis. It integrates multiple methods to fetch homologous templates, perform alignments, and calculate similarity scores. The tool supports flexible configurations for gap recognition, RMSD thresholds, and homology search methods.

---

## Features

- Aligns protein structures using user-specified or automated homology search methods.
- Flexible support for `.pdb` and `.cif` file formats.
- Gap recognition and RMSD filtering for precise alignments.
- Integration with PyMOL for visualization.
- Outputs results in a PyMOL session file (`.pse`) for further analysis.

---

## Installation

### 1. Clone the Repository

Clone the repository to your local machine:

```bash
git clone https://github.com/morth-lab/SIMalign.git
```

Navigate to the root of the repositoriy by:

```bash
cd SIMalign
```

### 2. Set Up the Conda Environment

Ensure you have `conda` installed, then use the provided `environment.yml` file to set up the environment:

```bash
conda env create -f environment.yml
conda activate simalign_env
```

<!-- ### 3. Install the Python Package

Install the SIMalign package using `pip`:

```bash
pip install .
```

This step will also download the `findSurfaceResidues.py` script into the `external/` directory. -->

---

## Usage

### Command-Line Interface

Run SIMalign using the following command:

```bash
python -m SIMalign.main --QUERY query.pdb
```

Where you substitute query.pdb with your input `.pdb` or `.cif` file.

### Available Arguments

| Argument                 | Description                                                                                                  |
|--------------------------|--------------------------------------------------------------------------------------------------------------|
| `--QUERY`                | Required. Path to the query protein structure file (`.pdb` or `.cif`).                                       |
| `--TEMPLATES`            | Optional. Folder with files used as templates for the query file (required for `user_specified` homology search method).             |
| `--HOMOLOGY_SEARCH_METHOD` | Method for homology search. Options: `foldseek`, `user_specified`, `BLASTp`. Default: `foldseek`.             |
| `--MAX_DISTANCE`         | Threshold distance for gap recognition. Default: `6`.                                                       |
| `--MAX_RMSD`             | Maximum allowed RMSD for alignments. Default: `5`.                                                          |
| `--FOLDSEEK_DATABASES`   | Databases for Foldseek. Options: `afdb50`, `afdb_swissprot`, `afdb_proteome`, `pdb100`. Default: `afdb50`.    |
| `--RESULT_DIR`           | Path to the directory for saving results. Default: `./<job_key>`.                                           |
| `--BLOSUM`               | Blosum

Folder with files used as templates for the query file (only required for user-specified homology search method).

---

## Output

- Results are saved in a directory named after the job key.
- Includes a PyMOL session file (`SIMalign_<job_key>.pse`) for further analysis.

---

## External Dependencies

This project requires the `findSurfaceResidues.py` script from the [PyMOL Script Repository](https://github.com/Pymol-Scripts/Pymol-script-repo).

The script is automatically downloaded during the `pip install .` step and stored in the `external/` directory.

---

## Contributing

We welcome contributions! Please fork the repository, make your changes, and submit a pull request.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- The PyMOL Script Repository for providing useful scripts like `findSurfaceResidues.py`.
- Code review and suggestions were enhanced using ChatGPT (January 2025).






<img src="logo.png" alt="Alt Text" width="200">

### [Publication](https://services.healthtech.dtu.dk/services/SIMAlign-1.0/)

### [Webserver](https://services.healthtech.dtu.dk/services/SIMAlign-1.0/)

<details><summary><b>Citation</b></summary>

If you use this code or the models in your research, please cite the following paper:

```bibtex
TEST TEST TEST
@article{rotilio2024structural,
  title={Structural and Functional Characterization of an Amidase Targeting a Polyurethane for Sustainable Recycling},
  author={Rotilio, Laura and Bayer, Thomas and Meinert, Hannes and Teixeira, Luis MC and Johansen, Martin B and Sommerfeldt, Andreas and Petersen, Allan R and Sandahl, Alexander and Keller, Malene B and Holck, Jesper and others},
  journal={Angewandte Chemie},
  pages={e202419535},
  year={2024},
  publisher={Wiley Online Library},
  doi={10.1002/anie.202419535}
}
```

</details>

<details open><summary><b>Table of contents</b></summary>

- [Usage](#usage)
  - [Setup Conda Environment](#environment)
  - [Run SIMalign](#SIMalign)
</details>


## Usage  <a name="usage"></a>

You can use the program directly on our [website](https://services.healthtech.dtu.dk/services/SIMAlign-1.0/) or run it locally using [Anaconda](https://docs.anaconda.com/anaconda/install/index.html).

### Setup Conda Environment <a name="environment"></a>

Clone the repository

	git clone https://github.com/morth-lab/SIMalign

Navigate to the root of the repository and run the following commands:

	conda env create --file environment.yml
	conda activate simalign

### Run SIMalign <a name="SIMalign"></a>

To run SIMalign, run the following command:

    python SIMalign.py -q <path_to_query> [optional_arguments]

It is required to input the path to the protein structure query file (.pdb or .cif). The entire list of options can be found using the help command:

    python SIMalign.py -h


To run SIMalign as a commandline write:
	$ python SIMALign_parser --QUERY path/to/query_file --TEMPLATE path/to/template_file_1 path/to/template_file_2 [...] --JOB_KEY job_key --HOMOLOGY_SEARCH_METHOD method --MAX_DISTANCE int --MAX_INITIAL_RMSD int --afdb50 True --afdb_swissprot True --afdb_proteome True --pdb100 True --FOLDSEEK_MODE "tmalign" --THRESHOLD float --NUMB_HOMO int --SEQUENCE_IDENTITY float --RESULT_DIR path/to/result_directory


int = an integer
float = a float number

The Following of the variables has default values: 
<variable>               == <default value>
--TEMPLATE               == None 
--JOB_KEY                ==  program-generated job-key 
--HOMOLOGY_SEARCH_METHOD == foldseek 
--MAX_DISTANCE           == 6 
--MAX_INITIAL_RMSD       == 5
--afdb50                 == True 
--afdb_swissprot         == False 
--afdb_proteome          == False 
--pdb100                 == False 
--FOLDSEEK_MODE          == tmalign 
--THRESHOLD              == 0.7 
--NUMB_HOMO              == 0 
--SEQUENCE_IDENTITY      == 0.6 
--RESULT_DIR             == User_data/{JOB_KEY}

the only obligatory option is --QUERY. 
When defining route to file, be sure to include the whole path from root. 


The program uses APIs (Foldseek ..) and therefore a prober internet connection is required for it to work optimal.

#################################
### Data Structure of Results ###
#################################

Then the program is done, the user should have a result directory of the following structure: 

RESULT_DIR:
|
|- SIMAlign_{JOB_KEY}.zip
|  |- tar.gz
|     |- m8-file
|  |- alignment.aln
|  |- m8-file
|  |- hotspots_mode_1.html
|  |- hotspots_mode_2.html
|  |- pymol_output.pse
|  |- AF-<x>-F1-model-v4.pdb
|  |-  <-||-> *x




## External Dependencies

The project depends on the `findSurfaceResidues.py` script from the [PyMOL Script Repository](https://github.com/Pymol-Scripts/Pymol-script-repo).

### Setup Instructions

1. During installation, the script will be automatically downloaded to the `external/` directory.
2. If needed, you can manually download it using:
   ```bash
   python setup_externals.py
