import argparse
import os
from core import SIMalign
import sys
from utils import validate_structure_file, encrypt_key, create_output_dirs, log_message, detect_structure_format
import shutil

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
        default=5,
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
    parser.add_argument(
        "--only_core",
        type=str,
        default="1",
        choices=["0", "1"],
        help="If set to 1, only hotspots in the core of the protein will be considered. If set to 0, all hotspots will be considered. Default is 1."
    )

    args = parser.parse_args()

    # Change "0" extension to ".pdb" for web server
    if args.QUERY.endswith(".0"):
        old_query_path = args.QUERY
        args.QUERY = detect_structure_format(old_query_path)
        # args.QUERY = args.QUERY[:-1] + "pdb"
        new_query_path = args.QUERY
        if new_query_path.endswith("0"):
            print(f'<p style="color:red;"><b>ERROR:</b> Could not detect structure format for query file ({old_query_path}). Please make sure the file has a valid structure format extension (.pdb or .cif).</p>')
            sys.exit(1)
        os.rename(old_query_path, new_query_path)


                

    # Validate arguments and inputs
    if not validate_structure_file(args.QUERY):
        print(f'<p style="color:red;"><b>ERROR:</b> Could not open or read query file</p>')
        sys.exit(1)


    if args.HOMOLOGY_SEARCH_METHOD == "user_specified":
        if args.TEMPLATES is None and args.TEMPLATES_DIR is None:
            print(f'<p style="color:red;"><b>ERROR:</b> Please provide either a list of template files or a directory containing template files when using user-specified homology search method.</p>')
            sys.exit(1)
        elif args.TEMPLATES is not None:
            templates = args.TEMPLATES
        else:
            templates = [os.path.join(args.TEMPLATES_DIR, temp_file) for temp_file in os.listdir(args.TEMPLATES_DIR)]
        if len(templates) < 2:
            print(f'<p style="color:red;"><b>ERROR:</b> Please provide 2 or more template files when using user-specified homology search method.</p>')
            sys.exit(1)

        # Change "0" extension to ".pdb" for web server
        for i, temp_file in enumerate(templates):
            print(temp_file)
            if temp_file.endswith(".0"):
                old_temp_file_path = temp_file
                new_temp_file_path = detect_structure_format(old_temp_file_path)
                temp_file = new_temp_file_path
                print(f"Detected template file format for {old_temp_file_path}: {new_temp_file_path}")
                if new_temp_file_path.endswith("0"):
                    templates.pop(i)
                    print(f'<p style="color:orange;"><b>WARNING:</b> Could not detect structure format for {old_temp_file_path}. This file will be skipped.</p>')
                    continue
                # new_temp_file_path = temp_file
                os.rename(old_temp_file_path, new_temp_file_path)
                templates[i] = new_temp_file_path
            if not validate_structure_file(temp_file):
                print(f'<p style="color:orange;"><b>WARNING:</b> Could not open or read {temp_file}</p>')
                templates.remove(temp_file)
        if len(templates) < 2:
            print(f'<p style="color:red;"><b>ERROR:</b> After validating the template files, less than 2 valid template files remain. Please provide at least 2 valid template files.</p>')
            sys.exit(1)
    else:
        templates = None


    if args.MAX_DISTANCE <= 0:
        print(f'<p style="color:red;"><b>ERROR:</b> The maximum sequence distance parameter (MAX_DISTANCE) must be greater than 0.</p>')
        sys.exit(1)

    if args.MAX_RMSD <= 0:
        print(f'<p style="color:red;"><b>ERROR:</b> The RMSD threshold (MAX_RMSD) must be greater than 0.</p>')
        sys.exit(1)

    if not (0.0 <= args.FOLDSEEK_THRESHOLD <= 1.0):
        print(f'<p style="color:red;"><b>ERROR:</b> The Foldseek threshold (FOLDSEEK_THRESHOLD) must be between 0.0 and 1.0.</p>')
        sys.exit(1)

    if args.NUMB_TEMPLATES < 2:
        print(f'<p style="color:red;"><b>ERROR:</b> Please use 2 or more template files. NUMB_TEMPLATES must be greater than 2.</p>')
        sys.exit(1)


    # Find MUSCLE binary
    muscle_path = shutil.which("muscle")
    if not muscle_path:
        print(f'<p style="color:red;"><b>ERROR:</b> MUSCLE binary not found on PATH</p>')
        sys.exit(1)

    tmp_dir, result_dir = create_output_dirs(args.RESULT_DIR, args.TMP_DIR)
    zip_file_path = os.path.join(result_dir, f"{job_key}_simalign")
    os.makedirs(zip_file_path, exist_ok=True)
    log_file_path = os.path.join(zip_file_path, f"{args.JOB_KEY}_log.txt")

    settings = [
        "SIMalign run settings:",
        f"Query = {args.QUERY}",
        f"job_key = {args.JOB_KEY}",
        f"result_dir = {args.RESULT_DIR}",
        f"tmp_dir = {args.TMP_DIR}",
        f"templates = {templates}",
        f"homology_search_method = {args.HOMOLOGY_SEARCH_METHOD}",
        f"max_dist = {args.MAX_DISTANCE}",
        f"max_rmsd = {args.MAX_RMSD}",
        f"foldseek_databases = {args.FOLDSEEK_DATABASES}",
        f"foldseek_mode = {args.FOLDSEEK_MODE}",
        f"foldseek_threshold = {args.FOLDSEEK_THRESHOLD}",
        f"numb_templates = {args.NUMB_TEMPLATES}",
        f"BLOSUM = {args.BLOSUM}",
        f"only_core = {args.only_core}",
        f"muscle_path = {muscle_path}",
        f"log_file_path = {log_file_path}"
    ]
    log_message(log_file_path, "\n".join(settings))



    # Run SIMalign
    SIMalign(query=args.QUERY,
             job_key=args.JOB_KEY,
             result_dir=zip_file_path,
             tmp_dir=tmp_dir,
             templates=templates,
             homology_search_method=args.HOMOLOGY_SEARCH_METHOD,
             max_dist=args.MAX_DISTANCE,
             max_rmsd=args.MAX_RMSD,
             foldseek_databases=args.FOLDSEEK_DATABASES,
             foldseek_mode=args.FOLDSEEK_MODE,
             foldseek_threshold=args.FOLDSEEK_THRESHOLD,
             numb_templates=args.NUMB_TEMPLATES,
             BLOSUM=args.BLOSUM,
             only_core=args.only_core,
             muscle_path=muscle_path,
             log_file_path=log_file_path
            )


    
    shutil.make_archive(zip_file_path, 'zip', zip_file_path)
    log_message(log_file_path, f"Results saved zip file: {zip_file_path}.zip")

    document_root = "/var/www/services"
    download_path = zip_file_path.replace(document_root, "")
    print('<img src="https://raw.githubusercontent.com/morth-lab/SIMalign/main/logo.png" width="200">')
    print(f'<a href="{download_path}.zip" download>Download here</a>')


if __name__ == "__main__":
    main()