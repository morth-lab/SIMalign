import argparse
import os
from core import SIMalign
import sys
from utils import validate_structure_file, encrypt_key, create_output_dirs

def main():

    # Generate a default job key
    job_key = encrypt_key()

    # Set up argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-q",
        "--QUERY",
        type=str,
        required=True,
        help="Required argument! Path to the query protein structure file (.pdb or .cif)."
    )
    parser.add_argument(
        "-t",
        "--TEMPLATES",
        nargs="+",
        default=None,
        help="List of paths to files used as templates for the query file (only required for user-specified homology search method)."
    )  
    parser.add_argument(
        "-t-dir",
        "--TEMPLATES_DIR",
        type=str,
        default=None,
        help="Folder with files used as templates for the query file (only required for user-specified homology search method)."
    ) 
    parser.add_argument(
        "-j",
        "--JOB_KEY",
        type=str,
        default=job_key,
        help="The job key (job name). The generated files will be saved under this name."
    )
    parser.add_argument(
        "-H",
        "--HOMOLOGY_SEARCH_METHOD",
        type=str,
        default="foldseek",
        choices=['foldseek', 'user_specified'],
        help="Method used to fetch homologue templates for usage in the program. Choose either foldseek or user_specified. (default is foldseek)."
    )
    parser.add_argument(
        "-d",
        "--MAX_DISTANCE",
        type=int,
        default=7,
        help="Threshold distance for gap recognition. If an amino acid has more than this distance to the query structure, it is recognized as a gap."
    )
    parser.add_argument(
        "-r",
        "--MAX_RMSD",
        type=int,
        default=5,
        help="Max allowed RMSD ('Root-Mean-Square Deviation'). Structures aligning with a higher RMSD than the given threshold will be removed from the algorithm."
    )
    parser.add_argument(
        "-fd",
        "--FOLDSEEK_DATABASES",
        nargs="+",
        choices=["afdb50", "afdb_swissprot", "afdb_proteome", "pdb100"],
        default=["afdb50"],
        help="Foldseek databases (only necessary for foldseek search method). Default: afdb50."
    )
    parser.add_argument(
        "-fm",
        "--FOLDSEEK_MODE",
        type=str,
        default="tmalign",
        choices=["tmalign", "3diaa"],
        help="Foldseek mode - Either 'TM-align' or '3Di/AA' (only necessary for foldseek search method). Default: tmalign."
    )
    parser.add_argument(
        "-ft",
        "--FOLDSEEK_THRESHOLD",
        type=float,
        default=0.7,
        help="Foldseek threshold (either TM or E-value threshold depending on Foldseek mode)."
    )
    parser.add_argument(
        "-nt",
        "--NUMB_TEMPLATES",
        type=int,
        default=20,
        help="Number of top performing Foldseek homologs based on (TM or E-value depending on Foldseek mode) that is used in the SIMalign algorithm."
    )
    parser.add_argument(
        "-sident",
        "--SEQUENCE_IDENTITY",
        type=float,
        default=0.6,
        help="Min. sequence identity (between 0-1) for BLASTp."
    ) 
    parser.add_argument(
        "-scov",
        "--SEQUENCE_COV",
        type=float,
        default=0.6,
        help="Min. sequence coverage (between 0-1) for BLASTp."
    ) 
    parser.add_argument(
        "-e",
        "--E_VALUE",
        type=float,
        default=0.001,
        help="E-value threshold for BLASTp."
    )
    parser.add_argument(
        "-rt",
        "--REDUNDANCY_THRESHOLD",
        type=float,
        default=0.9,
        help="Redundancy threshold for removing similar sequences from the BLASTp MSA."
    )
    parser.add_argument(
        "-R",
        "--RESULT_DIR",
        type=str,
        default=os.path.join(".", f"{job_key}"),
        help="The path to the result folder."
    )
    parser.add_argument(
        "-tmp",
        "--TMP_DIR",
        type=str,
        default=os.path.join(".", "tmp"),
        help="The path to the tmp folder."
    )   
    parser.add_argument(
        "-b",
        "--BLOSUM",
        type=str,
        default="BLOSUM62",
        choices=["BLOSUM50","BLOSUM62"],
        help="BLOSUM matrix used for sequence alignemnt."
    )
    args = parser.parse_args()

    # Validate arguments and inputs
    if not validate_structure_file(args.QUERY):
        print(f"ERROR: The query file {args.QUERY} must be of either .pdb or .cif format.")
        sys.exit(1)


    if args.HOMOLOGY_SEARCH_METHOD == "user_specified":
        if args.TEMPLATES is None and args.TEMPLATES_DIR is None:
            print("ERROR: Please provide either a list of template files or a directory containing template files when using user-specified homology search method.")
            sys.exit(1)
        elif args.TEMPLATES is not None:
            templates = args.TEMPLATES
        else:
            templates = [os.path.join(args.TEMPLATES_DIR, temp_file) for temp_file in os.listdir(args.TEMPLATES_DIR)]
        if len(templates) < 2:
            print("ERROR: Please provide 2 or more template files when using user-specified homology search method.")
            sys.exit(1)
        for temp_file in templates:
            if not validate_structure_file(temp_file):
                print(f"ERROR: {temp_file} must be of either .pdb or .cif format.")
                sys.exit(1)


    if args.MAX_DISTANCE <= 0:
        print("ERROR: The maximum sequence distance parameter (MAX_DISTANCE) must be greater than 0.")
        sys.exit(1)

    if args.MAX_RMSD <= 0:
        print("ERROR: The RMSD threshold (MAX_RMSD) must be greater than 0.")
        sys.exit(1)

    if not (0.0 <= args.FOLDSEEK_THRESHOLD <= 1.0):
        print("ERROR: The Foldseek threshold (FOLDSEEK_THRESHOLD) must be between 0.0 and 1.0.")
        sys.exit(1)

    if args.SEQUENCE_IDENTITY < 0.0 or args.SEQUENCE_IDENTITY > 1.0:
        print("ERROR: The sequence identity (SEQUENCE_IDENTITY) must be between 0.0 and 1.0.")
        sys.exit(1)

    if args.NUMB_TEMPLATES < 2:
        print("ERROR: Please use 2 or more template files. NUMB_TEMPLATES must be greater than 2.")
        sys.exit(1)




    # Print settings
    print("Running SIMalign with the following settings:")
    print(f"Query = {args.QUERY},",
          f"job_key = {args.JOB_KEY},",
          f"result_dir = {args.RESULT_DIR}",
          f"tmp_dir = {args.TMP_DIR}",
          f"templates = {templates},",
          f"homology_search_method = {args.HOMOLOGY_SEARCH_METHOD},",
          f"max_dist = {args.MAX_DISTANCE},",
          f"max_rmsd = {args.MAX_RMSD},",
          f"foldseek_databases = {args.FOLDSEEK_DATABASES},",
          f"foldseek_mode = {args.FOLDSEEK_MODE},",
          f"foldseek_threshold = {args.FOLDSEEK_THRESHOLD},",
          f"numb_templates = {args.NUMB_TEMPLATES},",
          f"sequence_identity = {args.SEQUENCE_IDENTITY},",
          f"sequence_cov = {args.SEQUENCE_COV},",
          f"redundancy_threshold = {args.REDUNDANCY_THRESHOLD},",
          f"BLOSUM = {args.BLOSUM}",
          sep="\n")

    tmp_dir, result_dir = create_output_dirs(args.RESULT_DIR, args.TMP_DIR)

    # Run SIMalign
    SIMalign(query=args.QUERY,
             job_key=args.JOB_KEY,
             result_dir=result_dir,
             tmp_dir=tmp_dir,
             templates=templates,
             homology_search_method=args.HOMOLOGY_SEARCH_METHOD,
             max_dist=args.MAX_DISTANCE,
             max_rmsd=args.MAX_RMSD,
             foldseek_databases=args.FOLDSEEK_DATABASES,
             foldseek_mode=args.FOLDSEEK_MODE,
             foldseek_threshold=args.FOLDSEEK_THRESHOLD,
             numb_templates=args.NUMB_TEMPLATES,
             sequence_identity=args.SEQUENCE_IDENTITY,
             sequence_cov=args.SEQUENCE_COV,
             redundancy_threshold=args.REDUNDANCY_THRESHOLD,
             BLOSUM=args.BLOSUM)

if __name__ == "__main__":
    main()
