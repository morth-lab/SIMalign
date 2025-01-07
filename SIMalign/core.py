import os
from Bio import AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from itertools import combinations
from statistics import median
import pymol2
import shutil


from .utils import foldseek_API_search, loading_structures_to_pymol, super_impose_structures, calculate_similarity_score, update_alignment_in_pymol, get_neighborAA, finding_hotspots, save_hotspot, format_pymol
from .models import StructureFile, Structure



def SIMalign(query, templates, job_key, homology_search_method, max_dist, max_rmsd, 
             foldseek_databases, foldseek_mode, foldseek_threshold, numb_templates, sequence_identity, 
             result_dir, tmp_dir, BLOSUM):
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
     

    elif homology_search_method == "BLASTp":
        pass

    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        stored = pymol.stored
        print("Loading structures into PyMOL...")
        structures = loading_structures_to_pymol(structure_files, query, cmd, stored)

        print("Superimposing structures...")
        structures = super_impose_structures(structures, max_rmsd, cmd, stored)

        print("Calculating similarity scores...")
        structures, align = calculate_similarity_score(structures, max_dist, cmd, BLOSUM)

        alignment_file_name = os.path.join(result_dir,"alignment.aln")
        print("Updating sequence alignment...")
        update_alignment_in_pymol(structures, align, cmd)
        cmd.save(alignment_file_name, selection="aln")
        
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
        result_zip = os.path.dirname(result_dir)
        shutil.make_archive(result_zip, 'zip', result_dir)

        print(f"Results saved at {result_zip}.zip")
