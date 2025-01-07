from random import randrange
from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import os
import sys
import subprocess
import requests
import time
from .models import StructureFile, Structure
from scipy.spatial import cKDTree
from Bio.Align import substitution_matrices
import numpy as np
from itertools import combinations
from collections import defaultdict
from statistics import median

## Main functions

def encrypt_key():
    """Generate a random 6-character job key."""
    key = ""
    library = "aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWyYzZ0123456789"
    for _ in range(6):
        key += library[randrange(len(library))]
    return key

def validate_structure_file(file_path):
    """Validate if a file is a valid PDB or CIF structure file."""
    parsers = [PDBParser(), MMCIFParser()]
    for parser in parsers:
        try:
            parser.get_structure('structure', file_path)
            return True
        except Exception:
            continue
    return False

def create_output_dirs(result_dir):
    """Create output directories for results and temporary files."""
    try:
        if not os.path.exists(result_dir):
            os.makedirs(os.path.join(result_dir, 'tmp'), exist_ok=True)
            os.makedirs(os.path.join(result_dir, 'results'), exist_ok=True)
    except OSError as e:
        print(f"ERROR: Could not create directories in {result_dir}: {e}")
        sys.exit(1)


## Core functions


def download_AF_structure(name,outfolder):
    name = name.split(" ")[0].split(".")[0]
    type = "pdb"
    url = "https://alphafold.ebi.ac.uk/files/"+name+"."+type
    print(f"\tDownloading {name}")
    os.system(f"wget -P {outfolder} {url}")
    return os.path.join(outfolder, name+"."+type)


def threeletter2oneletter(AA):
    """
    Convert three letter amino acid to one letter amino acid
    """
    try:
        amino_acid_translation = {
            'ALA': 'A',
            'ARG': 'R',
            'ASN': 'N',
            'ASP': 'D',
            'CYS': 'C',
            'GLU': 'E',
            'GLN': 'Q',
            'GLY': 'G',
            'HIS': 'H',
            'ILE': 'I',
            'LEU': 'L',
            'LYS': 'K',
            'MET': 'M',
            'PHE': 'F',
            'PRO': 'P',
            'SER': 'S',
            'THR': 'T',
            'TRP': 'W',
            'TYR': 'Y',
            'VAL': 'V',
            'MSE': 'M'}
        return amino_acid_translation[AA]
    except:
        return None

def aa_to_blosum_score(resn1,resn2,BLOSUM):
    index1 = BLOSUM.alphabet.index(threeletter2oneletter(resn1))
    index2 = BLOSUM.alphabet.index(threeletter2oneletter(resn2))
    return BLOSUM[index1,index2]

def average_coordinate(list_of_coordinates):
    """
    Get average coordinate of a list of coordinates
    """
    n = len(list_of_coordinates)
    for i, x in enumerate(list_of_coordinates):
        if i == 0:
            average = np.array(x)
        else:
            average += np.array(x)
    average = average/n
    return average


def process_nested_dicts(dicts, keys):
    """
    Optimized version of process_nested_dicts.
    """
    num_dicts = len(dicts)
    # Pre-filter dictionaries that lack required keys
    filtered_dicts = [
        (idx, d) for idx, d in enumerate(dicts)
        if all(key in d for key in keys)
    ]
    to_delete = defaultdict(set)

    # Compare filtered dictionaries pairwise
    for (idx_a, dict_a), (idx_b, dict_b) in combinations(filtered_dicts, 2):
        for key1, key2 in combinations(keys, 2):
            # Directly compare values (no redundant checks)
            if dict_a[key1] < dict_b[key1] and dict_a[key2] > dict_b[key2]:
                # Mark the appropriate dictionary for key deletion
                target_idx = idx_b if abs(dict_b[key1] - dict_b[key2]) > num_dicts else idx_a
                to_delete[target_idx].add(key2)

    # Remove keys in a single pass
    for idx, keys_to_remove in to_delete.items():
        dicts[idx] = {k: v for k, v in dicts[idx].items() if k not in keys_to_remove}

    return dicts

def dist_points(coord1, coord2):
    """
    Calculate distance between two points
    """
    coord1 = np.array(coord1)
    coord2 = np.array(coord2)
    distance = np.linalg.norm(coord2 - coord1)
    return distance


# Foldseek functions

def extract_highest_results(tresshold, number_of_templates, tmp_folder, result_folder, foldseek_mode):
    """
    Extract highest results from foldseek

    :param tresshold: tresshold for foldseek
    :param number_of_homologs: number of homologs to extract
    :param structure_files: list of structures to add foldseek results to
    :param foldseek_folder: folder where foldseek result is stored
    :param foldseek_mode: foldseek mode (tmalign or 3diaa)
    :return: list of structure_files with downloaded AF structures
    """
    for file in os.listdir(tmp_folder):
        m8file = os.path.join(tmp_folder,file)
    with open(m8file,"r") as infile:
        lines = infile.readlines()
        lines = lines[:number_of_templates]
        template_files = []
        for line in lines:
            line = line.split("\t")
            name = line[1].split(" ")[0]
            score = line[11]
            if (float(score) > tresshold and foldseek_mode == "tmalign") or (float(score) < tresshold and foldseek_mode == "3diaa"):
                file_location = None
                if name.startswith("AF"):
                    file_location = download_AF_structure(name,result_folder)
                    type = os.path.basename(file_location).split(".")[-1]
                else:
                    name = name.split("-")[0]
                    type = "PDB"
                template_files.append(StructureFile(name, file_location, type))
            elif len(template_files) < 2:
                print("\033[1;35mError: Not enough template structure files were downloaded!\033[0m")
                break
    return template_files


def foldseek_API_search(foldseek_mode, foldseek_databases, query, result_dir, foldseek_threshold, numb_templates):
    databases = ""
    for database in foldseek_databases:
        databases += f' -F "database[]={database}"'
    foldseekAPI = f'curl -X POST -F "q=@{query}" -F "mode={foldseek_mode}"{databases} https://search.foldseek.com/api/ticket'
    result = subprocess.run(foldseekAPI, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        ticket = result.stdout.split('"')[3]
    except:
        result = subprocess.run(foldseekAPI, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        ticket = result.stdout.split('"')[3]
    url = f"https://search.foldseek.com/api/ticket/{ticket}"
    # To check status of run
    response = requests.get(url)
    web_content = response.text
    # wait until foldseek search is complete:
    while not web_content.endswith('"COMPLETE"}\n'):
        time.sleep(2)
        response = requests.get(url)
        web_content = response.text

    # download foldseek result
    url = f"https://search.foldseek.com/api/result/download/{ticket}"
    tmp_folder = os.path.join(result_dir,"tmp")
    result_folder = os.path.join(result_dir,"results")
    os.system(f"wget {url} -O {tmp_folder}.tar.gz")
    os.system(f"tar -xzvf {tmp_folder}.tar.gz -C {tmp_folder}")
    os.system(f"rm {tmp_folder}.tar.gz")
    return extract_highest_results(float(foldseek_threshold), int(numb_templates), tmp_folder, result_folder, foldseek_mode)   

# PyMol functions

def loading_structures_to_pymol(structure_files,query,cmd,stored):
            
    cmd.reinitialize()
    for file in structure_files:
        # Extract attributes from the StructureFile object
        name = file.name
        file_location = file.file_location
        
        # Fetching PDB files
        if name == "PDB":
            try:
                cmd.fetch(name)
                print(f"\tFetched: {name}")
                Structure(name,cmd,stored).test_struc_format()
            except Exception as e:
                print(f"Fetch ERROR: {name} could not be fetched to PyMOL. Error: {e}")

        # Loading structures from paths
        else:
            # Normalize the file path
            if file_location not in query:
                normalized_path = os.path.normpath(file_location)
            else:
                normalized_path = os.path.normpath(query)
            
            # Check if the file exists at the specified path
            if not os.path.isfile(normalized_path):
                print(f"Error: File not found at {normalized_path}")
                continue
            
            try:
                # Load the structure in PyMOL using the normalized path and name
                cmd.load(normalized_path, name)
                print(f"\tLoaded: {name} from {normalized_path}")
                Structure(name,cmd,stored).test_struc_format()
            except Exception as e:
                print(f"Import ERROR: {normalized_path} could not be loaded in PyMOL. Error: {e}")

    structures = []
    for struc in cmd.get_object_list():
        structures.append(Structure(struc,cmd,stored))

    return structures


def super_impose_structures(structures, max_rmsd, cmd, stored):
    """
    Super impose structures to query_structure

    """
    
    n_templates = len(structures) - 1
    query_structure = structures[0]
    
    # Superimpose of all structures to query_structure
    for struc in structures[1:]:
        super = cmd.super(target=query_structure.first_chain,mobile=struc.first_chain)
        if super[0] > max_rmsd:
            print(f"\033[35m\tStructure {struc.name} was deleted because the RMSD to {query_structure.name} was above {max_rmsd}Å: {super[0]}Å\033[0m")
            cmd.delete(struc.name)
        elif super[1] < 100:
            print(f"\033[35m\tStructure {struc.name} was deleted because number of aligned atoms to {query_structure.name} was below 100: Aligned atoms={super[1]}\033[0m")
            cmd.delete(struc.name)
        else:
            print(f"\t{struc.name} was superimposed with: RMSD={round(super[0],3)}Å, Aligned atoms={super[1]}")

    # Update structures after super impostion
    structures = []
    for struc in cmd.get_object_list():
        structures.append(Structure(struc,cmd,stored))

    return structures



def calculate_similarity_score(structures, max_dist, cmd, BLOSUM_string):

    for struc in structures:
        struc.model = cmd.get_model(struc.CA)
        struc.cKDTree = cKDTree([atom.coord for atom in struc.model.atom])

    n_templates = len(structures) - 1
    
    BLOSUM = substitution_matrices.load(BLOSUM_string)
    align = []
    # LOOP through models
    for i, ref_struc in enumerate(structures):
        ref_kd = ref_struc.cKDTree

        # LOOP through CA atoms to find scores of one structure
        for ref_atom in ref_struc.model.atom:
            score = 0
            ref_resn = ref_atom.resn
            ref_coord = ref_atom.coord
            max_score = aa_to_blosum_score(ref_resn,ref_resn,BLOSUM) - np.min(BLOSUM)
            
            tmp_coordinates = [ref_coord]
            # Finding closest AA from other models to ref AA
            for struc in structures:
                if ref_struc != struc:
                    closest_pair = struc.cKDTree.query(ref_coord)
                    atom = struc.model.atom[closest_pair[1]]

                    # Validate that the atom pair is each others closest neighbor and below max_dist
                    if closest_pair[0] == ref_struc.cKDTree.query(atom.coord)[0] and closest_pair[0] <= max_dist:
                        score += ((aa_to_blosum_score(ref_resn,atom.resn,BLOSUM) - np.min(BLOSUM))/(n_templates*max_score))
                        tmp_coordinates.append(atom.coord)
            ref_struc.score_list.append(score)

            # Updating alignment based on structural infomation
            tmp_center = average_coordinate(tmp_coordinates)
            tmp = {struc.name: struc.cKDTree.query(tmp_center)[1]
                    for struc in structures if struc.cKDTree.query(tmp_center)[0] <= max_dist}

            
            if not any(all(tmp.get(struc.name) == ele.get(struc.name) for struc in structures if struc.name in ele) for ele in align):
                for x, ele in enumerate(align):
                    key = next(iter(tmp))  # Grab any key from tmp
                    if tmp.get(key) < ele.get(key, float('inf')):
                        align.insert(x, tmp)
                        break
                else:
                    align.append(tmp)


    return structures, align


def update_alignment_in_pymol(structures, align, cmd):
    """
    Update alignment in pymol
    """
    print("før",time.time())
    align = process_nested_dicts(align,[struc.name for struc in structures])
    print("efter",time.time())
    # Updating the align structure so it fits pymol
    new_align = []
    for pos in align:
        tmp = []
        for k,v in pos.items():
            tmp.append((k,int(v+1)))
        new_align.append(tmp)
    print(new_align)
    cmd.set_raw_alignment("aln",new_align)
    print(time.time())


# Functions for hotspot finding

def index_to_resi(index,align_seq,atomsCA):
    count = 0
    for i, AA in enumerate(align_seq.seq):
        if AA != "-":
            if i == index:
                if AA == "?":
                    return "-"
                else:
                    return int(atomsCA[count].resi)
            count += 1
        elif i == index:
            return "-"
        
def resi_to_index(residue,align_seq,atomsCA):
    count = 0
    for index, AA in enumerate(align_seq.seq):
        if AA != "-":
            if int(atomsCA[count].resi) == residue:
                return index
            count += 1

# def not_bigger_AA(ref_AA,target_AA):
#     """
#     Check if target AA is not bigger than ref AA
#     :param target_AA: target AA
#     :param ref_AA: ref AA
#     :return: True if target AA is bigger than ref AA, False otherwise
#     """
#     if ref_AA == target_AA:
#         return True
#     elif target_AA == "G":
#         return True
#     elif target_AA == "A" and ref_AA not in {"G","P"}:
#         return True
#     elif target_AA == "V" and ref_AA == "I":
#         return True
#     elif target_AA == "F" and ref_AA == "Y":
#         return True
#     elif target_AA == "L" and ref_AA in {"F","Y"}:
#         return True
#     elif target_AA == "S" and ref_AA in {"C","T"}:
#         return True
#     else:
#         return False
def bigger_AA(ref_AA,target_AA,if_refAA_and_targetAA_are_the_same=False):
    """
    Check if target AA is bigger than ref AA
    :return: True if target AA is bigger than ref AA, False otherwise
    :return: if_refAA_and_targetAA_is_the_same if target AA and ref AA are the same
    """
    if ref_AA == target_AA:
        return if_refAA_and_targetAA_are_the_same
    elif ref_AA == "G":
        return True
    elif ref_AA == "A" and target_AA not in {"G","P"}:
        return True
    elif ref_AA == "V" and target_AA == "I":
        return True
    elif ref_AA == "F" and target_AA == "Y":
        return True
    elif ref_AA == "L" and target_AA in {"F","Y"}:
        return True
    elif ref_AA == "S" and target_AA in {"C","T"}:
        return True
    else:
        return False

def median_40_percent(data):
    """
    Getting 40 percent median of the data
    """
    n = len(data)
    remove_count = int(0.2 * n)
    # Directly pass the sliced list to statistics.median()
    return median(sorted(data)[:n-remove_count])






def get_core(structures):
    """
    Get conserved core region
    :param structures: list of structures
    :return: list of AA in conserved region
    """
    score_list = structures[0].score_list
    med = median_40_percent(score_list)
    atoms = structures[0].model.atom
    # Get upper 60% median conserved AA
    for i, score in enumerate(score_list):
        if score > med:
            structures[0].conserved_AA.add(int(atoms[i].resi))

    # Add close AA within 6 AA to conserved_region
    conserved_region = structures[0].conserved_AA.copy()
    tmp = set()
    for atom in atoms:
        residue = int(atom.resi)
        if residue in conserved_region:
            tmp.update(set([int(atoms[x].resi) for x in structures[0].cKDTree.query_ball_point(atom.coord,6)]))

    conserved_region.update(tmp)

    # Getting SurfaceResidues:
    SurfaceResidues = set()
    external_dir = os.path.join(os.path.dirname(__file__), "..", "external")
    sys.path.append(external_dir)

    try:
        from findSurfaceResidues import findSurfaceResidues
        # Use functions or logic from the script
    except ImportError:
        raise ImportError("Could not import findSurfaceResidues. Ensure it is downloaded.")
        sys.exit(1)

    for residue in findSurfaceResidues(cmd, structures[0].first_chain):
        SurfaceResidues.add(residue[1])

    core = conserved_region - SurfaceResidues
    return list(core)

def update_neighborAA_list(neighborAA_list, tmp_neighborAAs, count, residue, resi_list, index, align_seq_query, atomsCA_query, align_seq_template, atomsCA_template):
    neighborAAs = set()

    # Adding glycines if any before because glycines don't have any side chain atoms
    while len(resi_list[:resi_list.index(residue)]) > count:
        neighborAA_list[index].append(set())  # Add empty set for glycines
        count += 1

    for ele in tmp_neighborAAs:
        neighborAAs.add(resi_to_index(ele, align_seq_query, atomsCA_query))  # Update close_AAs with the computed residue index
    neighborAAs.discard(resi_to_index(residue, align_seq_template, atomsCA_template))
    neighborAA_list[index].append(neighborAAs)  # Append the non-empty set
    count += 1  # Increase count only once after appending

    return neighborAA_list, count

def get_neighborAA(structures, align, cmd):
    """
    Get neighbor AA list and conserved regions
    :param structures: list of structures
    :param align: alignment
    :return: list of close AA
    """
    core = get_core(structures)
    query_atoms = cmd.get_model(structures[0].side_chains).atom
    query_kd = cKDTree([atom.coord for atom in query_atoms])

    core_index = [resi_to_index(residue, align[0], structures[0].model.atom) for residue in core]

    neighborAA_list = []
    for j, struc in enumerate(structures[1:]):
        resi_list = [index_to_resi(index, align[j+1], struc.model.atom) for index in core_index]
        resi_set = set(resi_list)

        neighborAA_list.append([])
        atoms = cmd.get_model(struc.side_chains).atom
        kd = cKDTree([atom.coord for atom in atoms])

        residue = 0
        tmp_close_AAs = set()
        count = 0

        for a in atoms:  # Remember, glycines have no side chains

            if residue != int(a.resi):
                if residue in resi_set:
                    neighborAA_list, count = update_neighborAA_list(neighborAA_list, tmp_close_AAs, count, residue, resi_list, j, align[0], structures[0].model.atom, align[j+1], struc.model.atom)
                    tmp_close_AAs = set()
                residue = int(a.resi)
            if residue in resi_set:
                # Finding atoms in query structure within 4Å of target template atom
                tmp_set = set([int(query_atoms[x].resi) for x in query_kd.query_ball_point(a.coord, 4)])
                tmp_close_AAs.update(tmp_set)

        if residue in resi_set:
            neighborAA_list, count = update_neighborAA_list(neighborAA_list, tmp_close_AAs, count, residue, resi_list, j, align[0], structures[0].model.atom, align[j+1], struc.model.atom)

        # Add empty set if sequence ends with glycine
        while len(resi_list) > count:
            neighborAA_list[j].append(set())
            count += 1

    return neighborAA_list, core, core_index



def finding_hotspots(neighborAA_list, align, structures, core, core_index, mode=1):
    """
    ChatGPT can you add some text here?
    """
    hotspot_list = []
    
    for query_atom in structures[0].model.atom:
        if int(query_atom.resi) in core:
            query_resi = int(query_atom.resi)
            query_resn = query_atom.resn
            ref_index = resi_to_index(query_resi, align[0], structures[0].model.atom)
            resi_list = []
            non_compatible_resi = []
            residue_to_mutate = set()
            resi_in_hotspot_list = set()
            
            # Loop through templates
            for j, seq in enumerate(align[1:]):
                resi_off = []
                k = j+1
                template_atoms = structures[k].model.atom # Template CA atoms
                template_AA = seq[ref_index] # Template AA
                non_compatible_tmpset = set()
                
                # Check that template AA is different than query AA and not gab
                if template_AA != "-" and template_AA != threeletter2oneletter(query_resn):
                    
                    # find atom that corresponds to template_AA
                    for template_atom in template_atoms:
                        if int(template_atom.resi) == index_to_resi(ref_index,seq,template_atoms):
                            break

                    # Check that CA atom of template is close to CA atom of query
                    if dist_points(query_atom.coord,template_atom.coord) < 1:    
                        # Check that ALL close AA to query AA is the same or smaller
                        for closeAA in neighborAA_list[j][core_index.index(ref_index)]:
                            if not bigger_AA(seq[closeAA],align[0][closeAA],if_refAA_and_targetAA_are_the_same=True):
                                resi_off.append(closeAA)
                            else:
                                non_compatible_tmpset.add(index_to_resi(closeAA,align[0],structures[0].model.atom))


                        # If searching for single mutations
                        if mode == 1:
                            # Check that ALL close AA to query AA is the same or smaller
                            if len(resi_off) == 0:
                                residue_to_mutate.add((query_resi, align[0][ref_index]))
                                # hotspot_flag = True
                                resi_list.append([(index_to_resi(ref_index,seq,template_atoms), seq[ref_index])])
                            else:
                                resi_list.append(["g"])

                        # If searching for double mutations
                        elif mode == 2:
                            # Check posiblity for double mutations
                            if len(resi_off) == 1 and query_resi not in resi_in_hotspot_list:
                                if resi_off[0] in core_index:
                                    resi_off_for_close = []
                                    for close_to_close in neighborAA_list[j][core_index.index(resi_off[0])]:
                                        if not bigger_AA(seq[close_to_close],align[0][close_to_close],if_refAA_and_targetAA_are_the_same=True):
                                            resi_off_for_close.append(close_to_close)
                                        else:
                                            non_compatible_tmpset.add(index_to_resi(close_to_close,align[0],structures[0].model.atom))
                                            
                                    if len(resi_off_for_close) == 1:
                                        if resi_off_for_close[0] == ref_index:
                                            residue_to_mutate.add((query_resi, align[0][ref_index]))
                                            residue_to_mutate.add((index_to_resi(resi_off[0],align[0],structures[0].model.atom), align[0][resi_off[0]]))
                                            resi_list.append([(index_to_resi(ref_index,seq,template_atoms), seq[ref_index]),  (index_to_resi(resi_off[0],seq,template_atoms), seq[resi_off[0]])])

                                        else:
                                            resi_list.append(["a"])
                                    elif len(resi_off_for_close) == 0:
                                        residue_to_mutate.add((query_resi, align[0][ref_index]))
                                        residue_to_mutate.add((index_to_resi(resi_off[0],align[0],structures[0].model.atom), align[0][resi_off[0]]))
                                        resi_list.append([(index_to_resi(ref_index,seq,template_atoms), seq[ref_index]),  (index_to_resi(resi_off[0],seq,template_atoms), seq[resi_off[0]])])
                                    else:
                                        resi_list.append(["b"])
                                else:
                                    resi_list.append(["c"])
                            else:
                                resi_list.append(["d"])
                    else:
                        resi_list.append(["e"])
                else:
                    resi_list.append(["f"])
                non_compatible_resi.append(non_compatible_tmpset)
                # Add hotspot to hotspot list
            if residue_to_mutate != set():
                hotspot_list.append([
                    residue_to_mutate,
                    resi_list,
                    non_compatible_resi
                    ])
                resi_in_hotspot_list.update(residue_to_mutate)
    return hotspot_list