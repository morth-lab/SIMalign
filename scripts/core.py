import os
from Bio import AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from itertools import combinations
from statistics import median
import pymol2
import shutil


from utils import foldseek_API_search, loading_structures_to_pymol, super_impose_structures, calculate_similarity_score, update_alignment_in_pymol, get_neighborAA, finding_hotspots, save_hotspot, format_pymol, save_scores_as_json

# Blast functions
#from utils import blastp, create_msa, remove_redundancy_from_msa, run_cd_hit

from utils import run_muscle, write_fasta, alignment_to_dict, update_msa

from models import StructureFile, Structure



def SIMalign(query, job_key, result_dir, tmp_dir="tmp", templates=None, homology_search_method="foldseek", max_dist=6, max_rmsd=5, 
             foldseek_databases=["afdb50"], foldseek_mode="tmalign", foldseek_threshold=0.7, numb_templates=20, sequence_identity=0.6, 
             sequence_cov=0.6, e_value=0.001, redundancy_threshold=0.9, BLOSUM="BLOSUM62"):
    """Run the SIMalign prediction algorithm."""

    basename = os.path.basename(query).split(".")
    query_file = StructureFile(basename[0], query, basename[1])
    structure_files = [query_file]

    if homology_search_method == "user_specified":
        for template in templates:
            basename = os.path.basename(template).split(".")
            structure_files.append(StructureFile(basename[0], template, basename[1]))

    elif homology_search_method == "foldseek":
        print("Running foldseek...")
        structure_files += foldseek_API_search(foldseek_mode, foldseek_databases, query, result_dir, tmp_dir, foldseek_threshold, numb_templates)
     

#    elif homology_search_method == "BLASTp":
#        print("Running BLASTp...")
#        with pymol2.PyMOL() as pymol:
#            cmd = pymol.cmd
#            stored = pymol.stored
#            query_sequence = loading_structures_to_pymol(structure_files, query, cmd, stored)[0].get_fasta()
#        blast_fasta = blastp(query_sequence, os.path.join(tmp_dir, "blastp.fasta"), database="nr", e_value=0.001, seq_identity=0.6, seq_coverage=0.6)
#        if blast_fasta:
#            non_redundant = run_cd_hit(blast_fasta, os.path.join(tmp_dir, "cd_hit.fasta"), threshold=redundancy_threshold)
#            if non_redundant:
#                ## BOLTZ-1
#                structure_files += [StructureFile(f"template_{i}", template, "fasta") for i, template in enumerate(non_redundant)]
#            # raw_msa_file = create_msa(blast_fasta, os.path.join(tmp_dir, "raw_msa.aln"))
#            # if raw_msa_file:
#            #     msa_file = remove_redundancy_from_msa(raw_msa_file, os.path.join(tmp_dir, "msa.aln"), threshold=redundancy_threshold)
#            #     if not msa_file:
#            #         raise RuntimeError("Could not filter MSA file.")
#            # else:
#            #     raise RuntimeError("Could not create MSA file.")
#            else:
#                raise RuntimeError("Could not remove redundancy from BLAST search.")
#        else:
#            raise RuntimeError("No hits were found from the BLAST search.")
#        ### BOLTZ-1


    
    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        stored = pymol.stored
        print("Loading structures into PyMOL...")
        structures = loading_structures_to_pymol(structure_files, query, cmd, stored)

        print("Superimposing structures...")
        structures = super_impose_structures(structures, max_rmsd, cmd, stored)

        print("Updating sequence alignment...")
        structure_names = []
        sequences = []
        for structure in structures:
            sequences.append(structure.get_fasta())
            structure_names.append(structure.name)
        sequences_path = os.path.join(tmp_dir, "sequences.fasta")
        write_fasta(sequences, structure_names, sequences_path)
        alignment_file_name = os.path.join(result_dir,"alignment.aln")
        run_muscle(sequences_path, alignment_file_name)
        align_dict = alignment_to_dict(alignment_file_name, structures, alignment_format="fasta")
        align_dict = update_msa(align_dict, structures, cmd, threshold=max_dist)



        update_alignment_in_pymol(align_dict, cmd, structure_names)
        cmd.save(alignment_file_name, selection="aln")

        print("Calculating similarity scores...")
        # Seqeunce alignment based on superimposition function from PyMOL
        structures, align = calculate_similarity_score(structures, max_dist, cmd, BLOSUM, alignment_file_name)



        save_scores_as_json(structures, result_dir)



        # alignment_file_name = os.path.join(result_dir,"alignment.aln")
        #print("Updating sequence alignment...")
        # Seqeunce alignment based on superimposition function from PyMOL
        # update_alignment_in_pymol(structures, align, cmd)
        
        
        print("Finding hotspots...")
        align = AlignIO.read(alignment_file_name,"clustal")
        neighborAA_list, core, core_index = get_neighborAA(structures, align, cmd)
            
        # Finding double mutations
        hotspot_list_double = finding_hotspots(neighborAA_list, align, structures, core, core_index, mode=2)
            
        # Finding single mutations
        hotspot_list_single = finding_hotspots(neighborAA_list, align, structures, core, core_index, mode=1)
            
        hotspot_list = hotspot_list_double + hotspot_list_single
        

        save_hotspot(hotspot_list_double, result_dir, structures, mode=2)
        save_hotspot(hotspot_list_single, result_dir, structures, mode=1)


        print("Coloring structures in PyMOL...")
        format_pymol(structures, hotspot_list, cmd)

        print("Saving results...")
        cmd.save(os.path.join(result_dir,"SIMalign_"+job_key+".pse"))
        # result_zip = os.path.dirname(result_dir)
        # shutil.make_archive(result_zip, 'zip', result_dir)

        print(f"Results saved in the folder: {result_dir}")
