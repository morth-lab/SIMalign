import os
from Bio import AlignIO
import pymol2


from utils import foldseek_API_search, loading_structures_to_pymol, super_impose_structures, calculate_similarity_score, update_alignment_in_pymol, get_neighborAA, finding_hotspots, save_hotspot, format_pymol, save_scores_as_json

from utils import run_muscle, write_fasta, alignment_to_dict, update_msa

from models import StructureFile, Structure



def SIMalign(query, job_key, result_dir, tmp_dir="tmp", templates=None, homology_search_method="foldseek", max_dist=6, max_rmsd=5, 
             foldseek_databases=["afdb50"], foldseek_mode="tmalign", foldseek_threshold=0.7, numb_templates=20, BLOSUM="BLOSUM62", only_core="1", muscle_path=None):
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
        sequences_path = os.path.join(result_dir, "sequences.fasta")
        write_fasta(sequences, structure_names, sequences_path)

        alignment_file_name = os.path.join(result_dir,"alignment.aln")

        print("Debug")
        print("cat")
        os.system(f"cat {sequences_path}")
        print("manual MUSCLE")
        os.system(f"/opt/conda/bin/muscle -align {sequences_path} -output {alignment_file_name}")

        run_muscle(sequences_path, alignment_file_name, muscle_path)
        print("Alignment file created:", alignment_file_name)
        align_dict = alignment_to_dict(alignment_file_name, structures, alignment_format="fasta")
        print("Updating MSA...")
        align_dict = update_msa(align_dict, structures, cmd, threshold=max_dist)
        print("Updating alignment in PyMOL...")

        update_alignment_in_pymol(align_dict, cmd, structure_names)
        cmd.save(alignment_file_name, selection="aln")

        print("Calculating similarity scores...")
        # Seqeunce alignment based on superimposition function from PyMOL
        structures, align = calculate_similarity_score(structures, max_dist, cmd, BLOSUM, alignment_file_name)

        save_scores_as_json(structures, result_dir)

        
        print("Finding hotspots...")
        align = AlignIO.read(alignment_file_name,"clustal")
        neighborAA_list, core, core_index = get_neighborAA(structures, align, cmd, only_core)
            
        # Finding double mutations
        hotspot_list_double = finding_hotspots(neighborAA_list, align, structures, core, core_index, only_core, mode=2)
            
        # Finding single mutations
        hotspot_list_single = finding_hotspots(neighborAA_list, align, structures, core, core_index, only_core, mode=1)
            
        hotspot_list = hotspot_list_double + hotspot_list_single
        

        save_hotspot(hotspot_list_double, result_dir, structures, mode=2)
        save_hotspot(hotspot_list_single, result_dir, structures, mode=1)


        print("Coloring structures in PyMOL...")
        format_pymol(structures, hotspot_list, cmd)

        print("Saving results...")
        cmd.save(os.path.join(result_dir,"SIMalign_"+job_key+".pse"))

        print(f"Results saved in the folder: {result_dir}")
