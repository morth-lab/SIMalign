# -*- coding: utf-8 -*-
from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import argparse
from random import randrange
import os
from SIMAlign_module_v23 import SIMAlign


def main():
    # stays as is. 
    def encrypt_key():
        
        key=""
        library="aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWyYzZ0123456789"
        for i in range(0,50):
            inte=randrange(0,len(library)-1)
            key+=library[inte]
        return key
    
    def is_integer(n): # yes, nitpick.
        try:
            float_n = float(n)
        except ValueError:
            return False
        else:
            return float_n.is_integer()

    # does this really need to be a function?
    # YES. It is needed for validation of more than one field. so yes. 

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
        help="the file that the algorithm should be run on (with path)")
    parser.add_argument(
        "-t",
        "--TEMPLATES",
        nargs="+",
        default=None,
        help="list of files with paths to files used as templates for query file in (only required for user-specified homology search method)"
        )    
    parser.add_argument(
        "-j",
        "--JOB_KEY",
        type=str,
        default=job_key,
        help="The job key that the files generated through this algorithm should be stored under."
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
        help="Max. Amino acid length for gap recognition" 
        )
    parser.add_argument(
        "-r",
        "--MAX_INITIAL_RMSD",
        type=int,
        default=5,
        help="Max. Allowed RMSD ('Root-Mean-Square deviation')" 
        ) 
    parser.add_argument(
        "--afdb50",
        type=bool,
        #action="store_true",
        choices=[True,False],
        help="search Foldseek using Alphafold/UniProt50 v4 (only necessary for foldseek search method)"
        ) 
    parser.add_argument(
        "--afdb_swissprot",
        type=bool,
        #action="store_false",
        choices=[True,False],
        help="search Foldseek using Alphafold/SwissProt v4 (only necessary for foldseek search method)"
        ) 
    parser.add_argument(
        "--afdb_proteome",
        type=bool,
        #action="store_false",
        choices=[True,False],
        help="search Foldseek using Alphafold/Proteome v4 (only necessary for foldseek search method)"
        ) 
    parser.add_argument(
        "--pdb100",
        type=bool,
        #action="store_false",
        choices=[True,False],
        help="search Foldseek using PDB100 (only necessary for foldseek search method)"
        )
    parser.add_argument(
        "--FOLDSEEK_MODE",
        type=str,
        default="tmalign",
        choices=["tmalign","3diaa"],
        help="Foldseek mode - Either 'TM-align' or '3Di/AA' (only necessary for foldseek search method)"
        )
    parser.add_argument(
        "--THRESHOLD",
        type=float,
        default=0.7,
        help="Threshold (either TM or E-value thresshold depending on Foldseek mode)"
        )
    parser.add_argument(
        "--NUMB_HOMO",
        type=int,
        default=20,
        help="Number of topperforming Foldseek homologs based on (TM or E-value depending on Foldseek mode)"
        )
    parser.add_argument(
        "--SEQUENCE_IDENTITY",
        type=float,
        default=0.6,
        help="Min. sequence identity (between 0-1) for BLASTp"
        ) 

    parser.add_argument(
        "-R",
        "--RESULT_DIR",
        type=str,
        default=f"./{job_key}",
        help="The working folder for the program."
        )
    args = parser.parse_args()

    ## checks for possible errors before running
    if not validate_structure_file(args.QUERY):
        print("ERROR:The Uploaded files must be of either .pdb or .cif format.")
        exit(1)

    if args.HOMOLOGY_SEARCH_METHOD=="user_specified":
        if args.TEMPLATES is None or len(args.TEMPLATES)<2: 
            print("ERROR:Please upload 2 or more files")
            exit(1)
        else:
            raise_error=False
            for i in args.TEMPLATES:
                if not validate_structure_file(i):
                    print(f"ERROR: {i} must be of either .pdb or .cif format.")
            if raise_error:
                print("ERROR: Uploaded template file(s) must be of either .pdb or .cif format.")
                exit(1)

    if args.MAX_DISTANCE<0:
        print("TThe maximum sequence distance parameter (MAX_DISTANCE) must be larger than 0.") 
        exit(1)                                                                
                                                                        
    if args.MAX_INITIAL_RMSD<0:
        print("The Value of MAX_INITIAL_RMSD cannot subceed 0.") 
        exit(1)                                                                   
    if args.THRESHOLD>1.0 or args.THRESHOLD<0:
        print("The Value of THRESHOLD must be between 0 and 1") 
        exit(1)
    if args.NUMB_HOMO<1:
        print("The number of homologues of (NUMB_HOMO) cannot subceed 1") 
        exit(1)                                                            
    if args.SEQUENCE_IDENTITY<0.0 or args.SEQUENCE_IDENTITY>1.0:
        print("minimum sequence identity (SEQUENCE_IDENTITY) must be between 0.0 and 1.0!") 
        exit(1)                                                    
    ## creates output-folder, if it does not exists.
    if not os.path.exists(args.RESULT_DIR):
        os.mkdir(args.RESULT_DIR)
    os.mkdir(args.RESULT_DIR+'/query')
    os.mkdir(args.RESULT_DIR+'/results')
    os.mkdir(args.RESULT_DIR+'/templates')
    
    
    try:
        print("query=",args.QUERY,
                "templates=",args.TEMPLATES,
                "job_key=",args.JOB_KEY,
                "homology_search_method=",args.HOMOLOGY_SEARCH_METHOD,
                "max_dist=",args.MAX_DISTANCE,
                "max_initial_rmsd=",args.MAX_INITIAL_RMSD,
                "afdb50=",args.afdb50,
                "afdb_swissprot=",args.afdb_swissprot,
                "afdb_proteome=",args.afdb_proteome,
                "pdb100=",args.pdb100,
                "foldseek_mode=",args.FOLDSEEK_MODE,
                "threshold=",args.THRESHOLD,
                "numb_Homo=",args.NUMB_HOMO,
                "sequence_identity=",args.SEQUENCE_IDENTITY,
                "result_dir=",args.RESULT_DIR)
        SIMAlign(query=args.QUERY,
                templates=args.TEMPLATES,
                job_key=args.JOB_KEY,
                homology_search_method=args.HOMOLOGY_SEARCH_METHOD,
                max_dist=args.MAX_DISTANCE,
                max_initial_rmsd=args.MAX_INITIAL_RMSD,
                afdb50=args.afdb50,
                afdb_swissprot=args.afdb_swissprot,
                afdb_proteome=args.afdb_proteome,
                pdb100=args.pdb100,
                foldseek_mode=args.FOLDSEEK_MODE,
                threshold=args.THRESHOLD,
                numb_Homo=args.NUMB_HOMO,
                sequence_identity=args.SEQUENCE_IDENTITY,
                result_dir=args.RESULT_DIR
        )
    except FileNotFoundError: 
        print(f"ERROR: {args.QUERY} not found")
    except Exception as e: # all other cases
        print("ERROR: an error occurred while performing the calculation:")
        print(str(e))
        print("Please contact the developers with input files and error message or open an issue on our GitHub repository https://www.github.com...")


if __name__ == "__main__":
    main()