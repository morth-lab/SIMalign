# -*- coding: utf-8 -*-
from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import argparse
from random import randrange
import os
from SIMalign_module import SIMalign
import sys


def main():

    def encrypt_key():
        
        key=""
        library="aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWyYzZ0123456789"
        for i in range(0,6):
            inte=randrange(0,len(library)-1)
            key+=library[inte]
        return key
    
    def is_integer(n): 
        try:
            float_n = float(n)
        except ValueError:
            return False
        else:
            return float_n.is_integer()

    def validate_structure_file(file_path):
        try:
            PDBParser().get_structure('structure',file_path)
            return True
        except:
            pass
        
        try:
            MMCIFParser().get_structure('structure',file_path)
            return True
        except:
            return False
    
    job_key=encrypt_key()
    
    parser=argparse.ArgumentParser()
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
    # parser.add_argument(
    #     "-t",
    #     "--TEMPLATE_DIR",
    #     type=str,
    #     default=None,
    #     help="Path to folder containing template files. Only used "
    #     )  
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
        choices=['foldseek','user_specified','BLASTp'],
        help="Method used to fetch homologue templates for usage in the program. Choose either foldseek, user_specified, or BLASTp. (default is foldseek)."
        )
    parser.add_argument(
        "-d",
        "--MAX_DISTANCE",
        type=int,
        default=6,
        help="Threshold distance for gap recognition. If an amino acid have more that this distance to the query structure it is recognized as a gap."
        )
    parser.add_argument(
        "-r",
        "--MAX_RMSD",
        type=int,
        default=5,
        help="Max allowed RMSD ('Root-Mean-Square deviation'). Structures aligning with a higher RMSD than the given threshold will be removed from the algorithm." 
        ) 
    parser.add_argument(
        "-fd",
        "--FOLDSEEK_DATABASES",
        nargs="+",
        choices=["afdb50","afdb_swissprot","afdb_proteome","pdb100"],
        default=["afdb50"],
        help="Foldseek databases (only necessary for foldseek search method). Default: afdb50."
        )
    parser.add_argument(
        "-fm",
        "--FOLDSEEK_MODE",
        type=str,
        default="tmalign",
        choices=["tmalign","3diaa"],
        help="Foldseek mode - Either 'TM-align' or '3Di/AA' (only necessary for foldseek search method). Default: tmalign."
        )
    parser.add_argument(
        "-ft",
        "--FOLDSEEK_THRESHOLD",
        type=float,
        default=0.7,
        help="Foldseek threshold (either TM or E-value thresshold depending on Foldseek mode)."
        )
    parser.add_argument(
        "--NUMB_HOMO",
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
        "-R",
        "--RESULT_DIR",
        type=str,
        default=os.path.join(".",f"{job_key}"),
        help="The path to the result folder."
        )
    args = parser.parse_args()

    ## checks for possible errors before running
    if not validate_structure_file(args.QUERY):
        print("ERROR: The query file must be of either .pdb or .cif format.")
        exit(1)

    if args.HOMOLOGY_SEARCH_METHOD=="user_specified":
        if args.TEMPLATES is None or len(args.TEMPLATES) < 2:
            print("ERROR: Please upload 2 or more template files when using user_specified homology search method.")
            exit(1)
        else:
            for temp_file in args.TEMPLATES:
                if not validate_structure_file(temp_file):
                    print(f"ERROR: {temp_file} must be of either .pdb or .cif format.")
                    exit(1)

    if args.MAX_DISTANCE<0:
        print("TThe maximum sequence distance parameter (MAX_DISTANCE) must be larger than 0.") 
        exit(1)                                                                                                                                 
    if args.MAX_RMSD<0:
        print("The Value of MAX_INITIAL_RMSD cannot subceed 0.") 
        exit(1)                                                                   
    if args.FOLDSEEK_THRESHOLD>1.0 or args.FOLDSEEK_THRESHOLD<0:
        print("The Value of THRESHOLD must be between 0 and 1.") 
        exit(1)
    if args.NUMB_HOMO<3:
        print("The number of homologues (NUMB_HOMO) should be at least 3.") 
        exit(1)                                                            
    if args.SEQUENCE_IDENTITY<0.0 or args.SEQUENCE_IDENTITY>1.0:
        print("minimum sequence identity (SEQUENCE_IDENTITY) must be between 0.0 and 1.0!") 
        exit(1)                     

    ## creates output-folder, if it does not exists.
    if not os.path.exists(args.RESULT_DIR):
        os.mkdir(args.RESULT_DIR)
    os.mkdir(os.path.join(args.RESULT_DIR,'query'))
    os.mkdir(os.path.join(args.RESULT_DIR,'results'))
    os.mkdir(os.path.join(args.RESULT_DIR,'templates'))
    
    # try:
    print("Running SIMalign with the following settings:")
    print(f"Query = {args.QUERY}",
            f"templates = {args.TEMPLATES}",
            f"job_key = {args.JOB_KEY}",
            f"homology_search_method = {args.HOMOLOGY_SEARCH_METHOD}",
            f"max_dist = {args.MAX_DISTANCE}",
            f"max_rmsd = {args.MAX_RMSD}",
            f"foldseek_databases = {args.FOLDSEEK_DATABASES}",
            f"foldseek_mode = {args.FOLDSEEK_MODE}",
            f"foldseek_threshold = {args.FOLDSEEK_THRESHOLD}",
            f"numb_homo = {args.NUMB_HOMO}",
            f"sequence_identity = {args.SEQUENCE_IDENTITY}",
            f"result_dir = {args.RESULT_DIR}",
            sep=", ")
    SIMalign(query=args.QUERY,
            templates=args.TEMPLATES,
            job_key=args.JOB_KEY,
            homology_search_method=args.HOMOLOGY_SEARCH_METHOD,
            max_dist=args.MAX_DISTANCE,
            max_initial_rmsd=args.MAX_RMSD,
            foldseek_databases=args.FOLDSEEK_DATABASES,
            foldseek_mode=args.FOLDSEEK_MODE,
            foldseek_threshold=args.FOLDSEEK_THRESHOLD,
            numb_Homo=args.NUMB_HOMO,
            sequence_identity=args.SEQUENCE_IDENTITY,
            result_dir=args.RESULT_DIR
    )
    # except FileNotFoundError: 
    #     print(f"ERROR: {args.QUERY} not found")
    # except Exception as e: # all other cases
    #     print("ERROR: an error occurred while performing the calculation:")
    #     print(str(e))
    #     print("Please contact the developers with input files and error message or open an issue on our GitHub repository https://www.github.com/morth-lab/SIMalign")


if __name__ == "__main__":
    main()