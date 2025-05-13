from random import randrange
from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import os
import sys
import subprocess
import requests
import time
import json
from models import StructureFile, Structure
from scipy.spatial import cKDTree
from Bio.Align import substitution_matrices
import numpy as np
from statistics import median
from Bio import AlignIO

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

def create_output_dirs(result_dir, tmp_dir):
    """Create output directories for results and temporary files."""
    try:
        if not os.path.exists(result_dir):
            os.makedirs(result_dir, exist_ok=True)
    except OSError as e:
        print(f"ERROR: Could not create directories in {result_dir}: {e}")
        sys.exit(1)
    try:
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir, exist_ok=True)
    except OSError as e:
        print(f"ERROR: Could not create directories in {tmp_dir}: {e}")
        sys.exit(1)
    return tmp_dir, result_dir


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




def dist_points(coord1, coord2):
    """
    Calculate distance between two points
    """
    coord1 = np.array(coord1)
    coord2 = np.array(coord2)
    distance = np.linalg.norm(coord2 - coord1)
    return distance


# Foldseek functions

def extract_highest_results(tresshold, number_of_templates, tmp_dir, result_dir, foldseek_mode):
    """
    Extract highest results from foldseek

    :param tresshold: tresshold for foldseek
    :param number_of_homologs: number of homologs to extract
    :param structure_files: list of structures to add foldseek results to
    :param foldseek_folder: folder where foldseek result is stored
    :param foldseek_mode: foldseek mode (tmalign or 3diaa)
    :return: list of structure_files with downloaded AF structures
    """
    for file in os.listdir(tmp_dir):
        if "_report" in file:
            continue
        m8file = os.path.join(tmp_dir,file)
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
                    file_location = download_AF_structure(name,result_dir)
                    type = os.path.basename(file_location).split(".")[-1]
                else:
                    name = name.split("-")[0]
                    type = "PDB"
                template_files.append(StructureFile(name, file_location, type))
            elif len(template_files) < 2:
                print("\033[1;35mError: Not enough template structure files were downloaded!\033[0m")
                break
    return template_files


def foldseek_API_search(foldseek_mode, foldseek_databases, query, result_dir, tmp_dir, foldseek_threshold, numb_templates):
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
    os.system(f"wget {url} -O {tmp_dir}.tar.gz")
    os.system(f"tar -xzvf {tmp_dir}.tar.gz -C {tmp_dir}")
    os.system(f"rm {tmp_dir}.tar.gz")
    return extract_highest_results(float(foldseek_threshold), int(numb_templates), tmp_dir, result_dir, foldseek_mode)   


# BLASTp functions


def run_muscle(input_fasta, output_fasta, muscle_cmd="muscle"):
    cmd = [
        muscle_cmd,
        "-align",   input_fasta,
        "-output",  output_fasta
    ]
    print("Running:", " ".join(cmd))
    # run it, raising an exception on error
    subprocess.run(cmd, check=True)


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
                Structure(name,cmd,stored).validate_structure_format()
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
                Structure(name,cmd,stored).validate_structure_format()
            except Exception as e:
                print(f"Import ERROR: {normalized_path} could not be loaded in PyMOL. Error: {e}")

    structures = []
    for struc in cmd.get_object_list():
        structures.append(Structure(struc,cmd,stored))

    return structures


def super_impose_structures(structures, max_rmsd, cmd, stored):
    """
    Superimpose structures to query_structure

    """
    
    n_templates = len(structures) - 1
    query_structure = structures[0]
    # cmd.alignto(query_structure.first_chain,object="aln")

    # Superimpose of all structures to query_structure
    for struc in structures[1:]:
        # super = cmd.super(target=query_structure.first_chain,mobile=struc.first_chain,object="aln")
        super = cmd.super(target=query_structure.first_chain, mobile=struc.first_chain)
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
        s = Structure(struc,cmd,stored)
        s.model = cmd.get_model(s.CA)
        s.cKDTree = cKDTree([atom.coord for atom in s.model.atom])
        structures.append(s)

    return structures



def calculate_similarity_score(structures, max_dist, cmd, BLOSUM_string, alignment_file_name):



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






def get_core(structures, cmd):
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
        from external.findSurfaceResidues import findSurfaceResidues

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
    core = get_core(structures, cmd)
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


def save_hotspot(hotspot_list, output_dir, structures, mode):
    doc_name=os.path.join(output_dir,"hotspots_mode_"+str(mode)+".html")
    doc=open(doc_name,"w")
    
    doc.write("------------------------------------------------------------------------------")
    if mode == 1:
        doc.write("<h1>Printing possible single mutations in "+structures[0].name+"</h1>")
    elif mode == 2:
        doc.write("<h1>Printing possible double mutations in "+structures[0].name+"</h1>")
    doc.write("------------------------------------------------------------------------------")
    doc.write("<h1>Legend</h1>")
    doc.write("<p style='color:rgb(209,14,14);'>Red : Mutations into smaller AAs</p>")
    doc.write("<p style='color:rgb(13,175,73);'>Green: Mutations into larger AAs</p>")
    doc.write("<p style='color:rgb(235,130,15);'>Orange: Undefined</p>")
    doc.write("------------------------------------------------------------------------------")
    printed = set()
    for hotspot in hotspot_list:

        # Finding unique mutations
        unique_mutations = {}
        
        for i, mutation_list in enumerate(hotspot[1]):
            if not isinstance(mutation_list[0], str):
                key = "|".join([x[1] for x in mutation_list])
                structure_name = structures[i+1].name

                non_comparable = set()
                for resi in hotspot[2][i]:
                    non_comparable.add(resi)

                # value = [structure_name, list(non_comparable)]
                if key in unique_mutations:
                    unique_mutations[key][0].append(structure_name)
                    unique_mutations[key][1].update(non_comparable)
                else:
                    unique_mutations[key] = [[structure_name], non_comparable]

        wt_list = hotspot[0]
        
        
        for k, v in unique_mutations.items():
            out_text = ""
            mutation_resn = k.split("|")
            mutated_structures = v[0]
            non_comparable = v[1]

            # color based on mutation size
            if all(bigger_AA(wt[1], mu) for wt, mu in zip(wt_list, mutation_resn)):
                out_text += "<p style='color:rgb(13,175,73);'>"
            elif all(bigger_AA(mu, wt[1]) for wt, mu in zip(wt_list, mutation_resn)):
                out_text += "<p style='color:rgb(209,14,14);'>"
            else:
                out_text += "<p style='color:rgb(235,130,15);'>"
            
        
            # write WT residues
            out_text += ", ".join([x[1]+str(x[0]) for x in wt_list])
            # write mutation residues
            out_text += " -> " + ", ".join(mutation_resn)
            # write non comparable residues
            if len(non_comparable) > 0:
                out_text += " | Non-comparable: " + ", ".join([str(x) for x in non_comparable])
            # write structure origins 
            out_text += " | Structures with given mutation: " + ", ".join(mutated_structures)
              
            doc.write(out_text)
    doc.write("</p>")
    doc.write("\n------------------------------------------------------------------------------")
    doc.close()

# pymol formatting functions

def color_by_number(number):
    """
    Gradient from red to white
    number 0: red
    number 1: white
    """
    if number >= 1:
        return [1,1,1]
    return [0.8+(number/5),number,number]

def color_structure(structure, cmd):
    """
    Coloring by similarity in pymol
    """
    score = structure.score_list
    cmd.set_color("unconserved", color_by_number(0))
    for s, atom in zip(score, structure.model.atom):
        if s == 0:
            cmd.color("unconserved", f"resi {atom.resi} and {structure.first_chain}")
        else:
            color_name = str(round(s,2))
            cmd.set_color(color_name, color_by_number(s))
            cmd.color(color_name, f"resi {atom.resi} and {structure.first_chain}")


def select_hotspots_in_pymol(hotspot_list, structures, cmd):
    count = 1
    for hotspot in hotspot_list:
        selection_string = f"({structures[0].first_chain} and ("
        resi_lists = [list(hotspot[0])]
        resi_lists.extend(hotspot[1])
        non_comparable = set()
        for non_comp in hotspot[2]:
            non_comparable.update(non_comp)
        if non_comparable != set():
            for non_comp in list(non_comparable):
                selection_string += f"resi {non_comp} or "
            selection_string = selection_string[:-4]+")) or ("
        else:
            selection_string = "("
        for i, resi_list in enumerate(resi_lists):
            if not isinstance(resi_list[0], str):
                selection_string += f"{structures[i].first_chain} and ("
                for resi in resi_list:
                    selection_string += f"resi {resi[0]} or "
                selection_string = selection_string[:-4]+")) or ("
        selection_string = selection_string[:-5]
        cmd.select(f"hotspot_{str(count)}", selection_string)
        count += 1
            



def format_pymol(structures, hotspot_list, cmd):
    cmd.hide("everything","all")
    for structure in structures:
        cmd.show("cartoon",structure.first_chain)
        print(f"\tColoring {structure.name}")
        color_structure(structure,cmd)
    cmd.hide("cgo", "aln")
    cmd.set("seq_view_label_mode", "1")
    select_hotspots_in_pymol(hotspot_list, structures, cmd)
    


def save_scores_as_json(structures, output_dir, filename="scores.json"):
    """
    Serialize a list of structures into a JSON file with this schema:
    {
      "structures": [
        {
          "name": <structure.name>,
          "residues": [<int>, <int>, ...],
          "scores": [<float>, <float>, ...]
        },
        ...
      ]
    }
    """

    # build up the list of structure dicts
    structs = []
    for struct in structures:
        # collect residue IDs and scores in parallel lists
        residues = [atom.resi for atom in struct.model.atom]
        scores   = list(struct.score_list)
        structs.append({
            "name":     struct.name,
            "residues": residues,
            "scores":   scores
        })

    # wrap it and write out
    output_path = os.path.join(output_dir, filename)
    with open(output_path, "w") as out:
        json.dump({"structures": structs}, out, indent=2)




def write_fasta(sequences: list, structure_names: list, filepath: str, line_width: int = 80) -> None:
    """
    Write a list of sequences to a FASTA file.

    Parameters
    ----------
    sequences : list
        sequence string (str).
    structure_names : list
        list of structure names (str).
    filepath : str
        Path to the output FASTA file.
    line_width : int, optional
        Maximum characters per line of sequence (default 80).
    """
    with open(filepath, "w") as out:
        for name, seq in zip(structure_names,sequences):
            out.write(f">{name}\n")
            # manually slice the sequence into chunks of line_width
            for i in range(0, len(seq), line_width):
                out.write(seq[i : i + line_width] + "\n")


def update_msa(msa_dicts, structures, cmd, threshold=6.0):
    """
    Refine an MSA (list of {struc_name: atom_index, ...} columns) by:
      1) Splitting any column whose members are >threshold Å apart into multiple sub-columns.
      2) Merging in any mutual nearest-neighbor pairs (within threshold) that aren't yet aligned.
    Returns a new list-of-dicts MSA.
    """


    index_to_pos_dicts = []
    for struc in structures:
        index_to_pos_dict = {}
        for i, atom in enumerate(struc.model.atom):
            index_to_pos_dict[atom.index] = i
        index_to_pos_dicts.append(index_to_pos_dict)


    # --- Split too-distant columns into clusters ---
    refined = []
    for col in msa_dicts:
        clusters = cluster_column(col, structures, threshold, index_to_pos_dicts)
        for cl in clusters:
            refined.append({s: r for s, r in cl})

    return refined



def distance(coord1, coord2):
    """Euclidean distance between two lists of coords."""
    return np.linalg.norm(np.array(coord1) - np.array(coord2))



def get_atom(structures, struc_name, atom_index, index_to_pos_dicts):
    for i, s in enumerate(structures):
        if s.name == struc_name:
            struc = s
            break
    else:
        raise ValueError(f"Structure {struc_name} not found.")
    
    # Get the atom using the index_to_pos_dict
    atom_pos = index_to_pos_dicts[i][atom_index]
    return struc.model.atom[atom_pos]



def cluster_column(col, structures, threshold, index_to_pos_dicts):
    """
    Cluster a single MSA column (dict struc->resi) so that every atom in each cluster
    is within `threshold` distance of every other (complete linkage).

    Args:
        col (dict): Mapping of structure names to residue indices.
        structures: Container of structure objects for atom lookup.
        threshold (float): Maximum allowed inter-atomic distance in a cluster.
        index_to_pos_dicts: Index-to-position lookup data for get_atom.

    Returns:
        List[List[(struc_name, atom_index)]]: Clusters of (structure, residue) pairs.
    """
    clusters = []
    centers = []  # numpy arrays of current cluster centroids

    for struc_name, atom_index in col.items():
        atom = get_atom(structures, struc_name, atom_index, index_to_pos_dicts)
        coord = atom.coord

        # Find nearest cluster center
        best_i = None
        best_dist = float('inf')
        for i, center in enumerate(centers):
            d = distance(coord, center)
            if d < best_dist:
                best_dist = d
                best_i = i

        # Try to add to the best cluster (complete linkage)
        if best_i is not None and best_dist <= threshold:
            candidate = clusters[best_i]
            if all(
                distance(coord, get_atom(structures, s, a, index_to_pos_dicts).coord) <= threshold
                for s, a in candidate
            ):
                candidate.append((struc_name, atom_index))
                # Update centroid
                all_coords = [get_atom(structures, s, a, index_to_pos_dicts).coord for s, a in candidate]
                centers[best_i] = np.mean(np.array(all_coords), axis=0)
                continue

        # Otherwise start a new cluster
        clusters.append([(struc_name, atom_index)])
        centers.append(np.array(coord))

    return clusters


def update_alignment_in_pymol(align_dict, cmd, keys):
    """
    Update alignment in pymol
    """

    # Updating the align structure so it fits pymol
    new_align = []
    for dict in align_dict:
        tmp = []
        for k in keys:
            if k in dict:
                tmp.append((k,dict[k]))
        new_align.append(tmp)
    cmd.set_raw_alignment("aln",new_align)



def alignment_to_dict(alignment_path, structures, alignment_format="fasta"):
    """
    Convert a multiple sequence alignment (MSA) to a list of dictionaries, where each dictionary
    represents a position in the MSA and contains the corresponding atom indices for each structure. 
    """
    
    alignment = AlignIO.read(alignment_path, alignment_format)

    structures_pos_added = {}
    for struct in structures:
        structures_pos_added[struct.name] = 0

    # build a mapping id → SeqRecord
    id2rec = { rec.description: rec for rec in alignment }

    msa_dicts = []
    for pos in range(len(alignment[0])):
        msa_dict = {}
        for struc in structures:
            seq = id2rec[struc.name].seq
            if seq[pos] != "-":
                atom_pos = structures_pos_added[struc.name]
                msa_dict[struc.name] = struc.model.atom[atom_pos].index
                structures_pos_added[struc.name] += 1
        msa_dicts.append(msa_dict)
    return msa_dicts