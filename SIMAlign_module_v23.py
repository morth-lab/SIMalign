import subprocess
import os
os.environ['PYMOL_GIT_MOD']='headless'
import subprocess
import time
import requests
from scipy.spatial import cKDTree
import numpy as np
from itertools import combinations
from Bio import AlignIO, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import statistics
import pymol2
import shutil
import sys

# the SIMAlign algorithm
def SIMAlign(query,templates,job_key ,homology_search_method,max_dist,max_initial_rmsd,afdb50,afdb_swissprot,afdb_proteome,pdb100,foldseek_mode,threshold,numb_Homo,sequence_identity,result_dir):
    # SIMAlign implementation with subprocess calls to PyMOL if necessary
    result_root=result_dir+"/"
    result_dir+="/results"
    
    foldseek_databases = []
    if afdb50:
      foldseek_databases.append("afdb50")
    if afdb_swissprot:
      foldseek_databases.append("afdb-swissprot")
    if afdb_proteome:
      foldseek_databases.append("afdb-proteome")
    if pdb100:
      foldseek_databases.append("pdb100")
    
    """# RUN prediction"""
    #@title Install
    #@markdown This may take a while.
    class StructureFile:
      def __init__(self, name, file_location, type):
        self.name = name
        self.file_location = file_location
        self.type = type

    def loading_file(path_list): # added by Alexander 18/10-24
    	if type(path_list)!=list:
    		path_list=[path_list]
    	files={}
    	for path in path_list:
    		name=path.split("/")[-1]
    		content=""
    		parsing_file=open(path,'r')
    		for line in parsing_file:
    			content+=line
    		files[name]=bytes(content,"utf-8")
    	return files
    
    query_file = list(loading_file(query))[0]
    query_file_list = [StructureFile(query_file.split(".")[0], query_file, query_file.split(".")[1])]
    if homology_search_method == "user_specified":
    	structure_files = query_file_list.copy()
    	for fn in templates:
            name=fn.split("/")[-1].split(".")[0]
            TYPE=fn.split(".")[-1]
            structure_files.append(StructureFile(name, fn, TYPE))
    
    def select_by_list(selection_list):
        """
        Make pymol selection by list of residues
        :param selection_list: list of residues
        :return: selection
        """
        out = " and ("
        for i, resi in enumerate(selection_list):
            out += f"resi {resi} or "
        return out[:-4]+")"
    
    #@title RUN SIMalign
   
    # -------------------------------- Foldseek -----------------------------------
    #def init_pymol():
    #    pymol.pymol_argv = ['pymol', '-cq']  # '-cq' to start without GUI and in quiet mode
    #    pymol.invoke()
    # ---------------- Foldseek functions -----------------
    #init_pymol()
    def make_empty_dir(folder_name):    
      """
      Make empty directory
      :param folder_name: name of the directory
      :return: None
      """
      if os.path.exists(folder_name):
        for file in os.listdir(folder_name):
          file_path = os.path.join(folder_name, file)
          os.system(f"rm {file_path}")
      else:
        os.system(f"mkdir {folder_name}")

    class Structure:
         def __init__(self, name, cmd):
             self.name = name
             self.cmd = cmd
             self.not_HETATM = self.name + " and not HETATM"
             self.first_chain = self.not_HETATM + " and chain " + cmd.get_model(self.not_HETATM).atom[0].chain
             self.CA = self.first_chain + " and name CA"
             self.side_chains = self.first_chain + " and (not name CA and not name N and not name C and not name O)"
             self.model = None
             self.cKDTree = None
             self.score_list = []
             self.conserved_AA = set()

         def get_canonicalAA_len(self):
            return len(cmd.get_model(self.CA).atom)

         def get_fasta_len(self):
             """
             Get fasta sequence from structure
             :param structure: structure to get fasta sequence from
             :return: fasta sequence
             """
             lines = cmd.get_fastastr(self.name).split("\n")
             fasta = ""
             for line in lines[1:]:
                 if line.startswith(">"):
                     return fasta
                 else:
                     fasta += line
             return len(fasta)

         def test_pdb_format(self):
             """
             Checking if the pdb format is valid
             """
             stored.residues = []
             cmd.iterate(self.name,"stored.residues.append(resi)")
             if [resi for resi in stored.residues if not resi.isdigit()]:
                 print("\033[35m\t"+self.name,"was removed due to error in its residues\033[0m")
                 cmd.delete(self.name)
             if self.get_fasta_len() != self.get_canonicalAA_len():
                 print("\033[35m\t"+self.name,"was removed because it contain non canonical amino acids\033[0m")
                 cmd.delete(self.name)
                 
    def extract_highest_results(tresshold, number_of_homologs, structure_files, foldseek_folder, foldseek_mode):
      """
      Extract highest results from foldseek

      :param tresshold: tresshold for foldseek
      :param number_of_homologs: number of homologs to extract
      :param structure_files: list of structures to add foldseek results to
      :param foldseek_folder: folder where foldseek result is stored
      :param foldseek_mode: foldseek mode (tmalign or 3diaa)
      :return: list of structure_files with downloaded AF structures
      """
      for file in os.listdir(foldseek_folder):
        m8file = os.path.join(foldseek_folder,file)
      with open(m8file,"r") as infile:
        lines = infile.readlines()
        lines = lines[:number_of_homologs]
        for line in lines:
          line = line.split("\t")
          name = line[1]
          score = line[11]
          if (float(score) > tresshold and foldseek_mode == "tmalign") or (float(score) < tresshold and foldseek_mode == "3diaa"):
            file_location = None
            if name.startswith("AF"):
              name = name.split(" ")[0].split(".")[0]
              type = "pdb"
              url = "https://alphafold.ebi.ac.uk/files/"+name+"."+type
              print(f"\tDownloading {name}")
              os.system(f"wget -P {foldseek_folder} {url}")
              file_location = os.path.join(foldseek_folder, name+"."+type)
            else:
              name = name.split("-")[0]
              type = "PDB"
            structure_files.append(StructureFile(name, file_location, type))
          elif len(structure_files) < 3:
            print("\033[1;35mError: Not enough structure files were downloaded!\033[0m")
            break
      return structure_files
    
    def aa_to_blosom(aa1,aa2):
        """
        Get blosom62 score of two amino acids
        :param aa1: first amino acid (Upper case)
        :param aa2: second amino acid (Upper case)
        :return: blosom62 score of two amino acids
        
        table from https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        """
        blosom62 = [
            [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4],
            [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4],
            [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4],
            [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4],
            [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
            [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4],
            [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
            [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4],
            [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4],
            [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4],
            [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4],
            [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4],
            [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4],
            [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4],
            [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
            [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4],
            [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4],
            [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4],
            [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4],
            [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4],
            [-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4],
            [-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
            [0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4],
            [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1]]
        aa_id_list = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","B","Z","X","-"]
        aa_to_id = dict()
        for i, aa in enumerate(aa_id_list):
            aa_to_id[aa] = i
        id1 = aa_to_id[aa1]
        id2 = aa_to_id[aa2]
        return blosom62[id1][id2]

    def average_coordinate(list_of_coordinates):
        """
        Get average coordinate of a list of coordinates
        :param list_of_coordinates: list of coordinates
        :return: average coordinate
        """
        n = len(list_of_coordinates)
        for i, x in enumerate(list_of_coordinates):
            if i == 0:
                average = np.array(x)
            else:
                average += np.array(x)
            average = average/n
        return average
    
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
                
    def threeletter2oneletter(AA):
        """
        Convert three letter amino acid to one letter amino acid
        :param AA: three letter amino acid
        :return: one letter amino acid
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
        
    def dist_points(coord1, coord2):
        """
        Calculate distance between two points
        :param coord1: first point
        :param coord2: second point
        :return: distance between two points
        """
        coord1 = np.array(coord1)
        coord2 = np.array(coord2)
        distance = np.linalg.norm(coord2 - coord1)
        return distance
    # Loop through core AAs in query structure
    def not_bigger_AA(ref_AA,target_AA):
        """
        Check if target AA is not bigger than ref AA
        :param target_AA: target AA
        :param ref_AA: ref AA
        :return: True if target AA is bigger than ref AA, False otherwise
        """
        if ref_AA == target_AA:
            return True
        elif target_AA == "G":
            return True
        elif target_AA == "A" and ref_AA not in {"G","P"}:
            return True
        elif target_AA == "V" and ref_AA == "I":
            return True
        elif target_AA == "F" and ref_AA == "Y":
            return True
        elif target_AA == "L" and ref_AA in {"F","Y"}:
            return True
        elif target_AA == "S" and ref_AA in {"C","T"}:
            return True
        else:
            return False
    
    def process_nested_dicts(dicts, keys):
        """
        Processes a list of dictionaries by comparing and potentially deleting key-value pairs
        based on the values of two keys from each dictionary.
        
        param: dicts (list): A list of dictionaries to process.
        param: keys (list): A list of two keys to compare and potentially delete.
        
        Returns: list: The processed list of dictionaries.
        """
        num_dicts = len(dicts)
        
        # Use indices instead of dictionaries as keys for the to_delete dictionary
        to_delete = {i: set() for i in range(num_dicts)}
        
        # Compare pairs of dictionaries
        for idx_a, idx_b in combinations(range(num_dicts), 2):
            dict_a, dict_b = dicts[idx_a], dicts[idx_b]
            for key1, key2 in combinations(keys, 2):
                if all(k in dict_a and k in dict_b for k in (key1, key2)):
                    if dict_a[key1] < dict_b[key1] and dict_a[key2] > dict_b[key2]:
                        # Use index instead of dictionary as the key in the to_delete dictionary
                        to_delete[idx_b if abs(dict_b[key1] - dict_b[key2]) > num_dicts else idx_a].add(key2)
                        
        # Remove the keys that were marked for deletion
        for idx, keys_to_remove in to_delete.items():
            for key in keys_to_remove:
                if key in dicts[idx]:
                    del dicts[idx][key]

        return dicts

    # -----------------------------------------------------
    with pymol2.PyMOL() as pymol:
        cmd=pymol.cmd
        stored=pymol.stored
        print("pymol session running for:",query)
        
        if homology_search_method == "foldseek":
            print("Running foldseek...")
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
            foldseek_folder=result_dir
            os.system(f"wget {url} -O {foldseek_folder}.tar.gz")
            make_empty_dir(result_dir)
            os.system(f"tar -xzvf {foldseek_folder}.tar.gz -C {foldseek_folder}")
            os.system(f"rm {foldseek_folder}.tar.gz")
            structure_files = query_file_list.copy()
            structure_files = extract_highest_results(float(threshold), int(numb_Homo), structure_files, result_dir, foldseek_mode)
     
        if homology_search_method == "BLASTp":
            print("Running BLASTp...")

            query_type=query.split(".")[1]
            file_id=query.split("/")[-1].split(".")[0]
            three_to_one = {
                'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
                'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                }
            fasta=""
            if query_type == "pdb":
                # Parse the PDB file
                with open(query, 'r') as file:
                    last_res=0
                    for line in file:
                        if line.startswith("ATOM"):
                            residue_number = int(line.split()[5])
                            residue_name = line.split()[3]
                            if residue_name in three_to_one and residue_number>last_res:
                                fasta += three_to_one[residue_name]
                                last_res=residue_number

            elif query_type == "cif":
                # Parse the CIF file
                with open(query, 'r') as file:
                    print("cif-file opened")
                    for line in file:
                        if line.startswith("_atom_site.label_comp_id"):
                            residue_name = line.split()[1]
                            if residue_name in three_to_one:
                                fasta += three_to_one[residue_name]
            print(fasta)
            #BLASTp_folder = "User_data/"+job_key+"/results"
            make_empty_dir(result_dir)
            structure_files = query_file_list.copy()
            def run_blast(sequence, key, identity_threshold=0.6, max_results=10):
                """
                Perform BLAST search and fetch homologous PDB structures.
                
                :param sequence: Amino acid sequence to BLAST.
                :param identity_threshold: Minimum sequence identity threshold.
                :param max_results: Maximum number of homologues to fetch.
                :return: List of PDB IDs matching the criteria.
                """
                # Step 2: Perform BLAST search
                print("Running BLAST search...")
                sequence=f">{key}\n{sequence}"
                print(sequence)
                result_handle = NCBIWWW.qblast("blastp", "pdb", sequence)

                # Step 3: Parse BLAST results
                print("Parsing BLAST results...")
                blast_records = NCBIXML.parse(result_handle)
                pdb_ids = []
                print(blast_records)
                for blast_record in blast_records:
                    print("processing ",blast_record)
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            identity = hsp.identities / hsp.align_length  # Calculate sequence identity
                            if identity >= identity_threshold:
                                pdb_id = alignment.accession
                                pdb_ids.append(pdb_id)
                                if len(pdb_ids) >= max_results:
                                    break
                        if len(pdb_ids) >= max_results:
                            break
                    if len(pdb_ids) >= max_results:
                        break
    
                print(f"Found {len(pdb_ids)} homologues with identity >= {identity_threshold}")
                return pdb_ids

            def fetch_pdb_files(pdb_ids, output_dir="pdb_files"):
                """
                Fetch PDB/CIF files for the given PDB IDs from RCSB PDB.
                
                :param pdb_ids: List of PDB IDs to fetch.
                :param output_dir: Directory to save the PDB files.
                """
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                    
                print(f"Fetching {len(pdb_ids)} PDB files...")
                base_url = "https://files.rcsb.org/download/"
                for pdb_id in pdb_ids:
                    pdb_url = f"{base_url}{pdb_id}.pdb"
                    response = requests.get(pdb_url)
                    if response.status_code == 200:
                        file_path = os.path.join(output_dir, f"{pdb_id}.pdb")
                        with open(file_path, "w") as f:
                            f.write(response.text)
                        print(f"Downloaded: {pdb_id}.pdb")
                        structure_files.append(StructureFile(pdb_id, output_dir, "pdb"))
                    else:
                        print(f"Failed to fetch: {pdb_id}.pdb")
   
            print("here we go !")
            homologues=run_blast(fasta, key=file_id, identity_threshold=sequence_identity)
            fetch_pdb_files(homologues,output_dir=result_dir)
            downloaded_files = os.listdir(result_dir)
            if not downloaded_files:
                raise FileNotFoundError("No PDB files were downloaded!")
            
            print(f"# Downloaded files: {downloaded_files}")

        # -------------------------------- SIMalign -----------------------------------

        # ---------------- SIMalign functions -----------------
        
        def loading_structures_to_pymol(structure_files,query,cmd):
                """
                Load structures to pymol

                :param structure_files: list of structures to load
                :return: query structure and list of structures
                """
                 
                cmd.reinitialize()
                for structure in structure_files:
                    # Extract attributes from the StructureFile object
                    name = structure.name
                    file_location = structure.file_location
                    
                    # Normalize the file path
                    if file_location not in query:
                        normalized_path = os.path.normpath(file_location)
                    else:
                        normalized_path=os.path.normpath(query)
                    
                    print("Loading structures to PyMOL...")
                    print(f"name: {name} ({type(name)})")
                    print(f"path: {normalized_path} ({type(normalized_path)})")
                    
                    # Check if the file exists at the specified path
                    if not os.path.isfile(normalized_path):
                        print(f"Error: File not found at {normalized_path}")
                        continue
                    
                    try:
                        # Load the structure in PyMOL using the normalized path and name
                        cmd.load(normalized_path, name)
                        print(f"Loaded: {name} from {normalized_path}")
                        Structure(name,cmd).test_pdb_format()
                    except Exception as e:
                        print(f"Import ERROR: {normalized_path} could not be loaded in PyMOL. Error: {e}")
   
                structures = []
                for struc in cmd.get_object_list():
                    structures.append(Structure(struc,cmd))

                return structures

        def super_impose_structures(structures, max_initial_rmsd,cmd):
                """
                Super impose structures to query_structure
                
                :param structures: list of structures to super impose
                :return: list of structures after super imposed
                """
                
                n_homologous_list = len(structures) - 1
                query_structure = structures[0]
                
                # Superimpose of all structures to query_structure
                for struc in structures[1:]:
                    super = cmd.super(target=query_structure.first_chain,mobile=struc.first_chain)
                    if super[0] > max_initial_rmsd:
                        print(f"\033[35m\tStructure {struc.name} was deleted because the RMSD to {query_structure.name} was above {max_initial_rmsd}Å: {super[0]}Å\033[0m")
                        cmd.delete(struc.name)
                    elif super[1] < 100:
                        print(f"\033[35m\tStructure {struc.name} was deleted because number of aligned atoms to {query_structure.name} was below 100: Aligned atoms={super[1]}\033[0m")
                        cmd.delete(struc.name)
                    else:
                        print(f"\t{struc.name} was superimposed with: RMSD={round(super[0],3)}Å, Aligned atoms={super[1]}")

                # Update structures after super impostion
                structures = []
                for struc in cmd.get_object_list():
                    structures.append(Structure(struc,cmd))

                return structures

        def SIMalign_algorithm(structures, max_dist):
                """
                SIMalign algorithm
                
                :param structures: list of structures to align
                :param max_dist: maximum distance between two amino acids before it is considered as a gab in the alignment
                :return: list of structures after alignment and scoring
                """

                for struc in structures:
                    struc.model = cmd.get_model(struc.CA)
                    struc.cKDTree = cKDTree([atom.coord for atom in struc.model.atom])

                n_homologous_list = len(structures) - 1

                align = []
                # LOOP through models
                for i, ref_struc in enumerate(structures):
                    ref_kd = ref_struc.cKDTree

                    # LOOP through CA atoms to find scores of one model
                    for ref_atom in ref_struc.model.atom:
                        score = 0
                        ref_resn = ref_atom.resn
                        ref_coord = ref_atom.coord
                        max_score = aa_to_blosom(ref_resn,ref_resn) + 4
                        
                        tmp_coordinates = [ref_coord]
                        # Finding closest AA from other models to ref AA
                        for struc in structures:
                            if ref_struc != struc:
                                closest_pair = struc.cKDTree.query(ref_coord)
                                atom = struc.model.atom[closest_pair[1]]
                                if closest_pair[0] == ref_struc.cKDTree.query(atom.coord)[0] and closest_pair[0] <= max_dist:
                                    score += ((aa_to_blosom(ref_resn,atom.resn) + 4)/(n_homologous_list*max_score))
                                    tmp_coordinates.append(atom.coord)
                        ref_struc.score_list.append(score)
                        tmp_center = average_coordinate(tmp_coordinates)
                        tmp = {struc.name: struc.model.atom[struc.cKDTree.query(tmp_center)[1]].index
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
            
        def update_alignment_in_pymol(structures, align):
                """
                Update alignment in pymol
                
                :param structures: list of structures to update alignment for
                :param align: list of alignment to update
                :return: None
                """

                align = process_nested_dicts(align,[struc.name for struc in structures])
                new_align = []
                for pos in align:
                    tmp = []
                    for k,v in pos.items():
                        tmp.append((k,v))
                    new_align.append(tmp)
                cmd.set_raw_alignment("aln",new_align)

            # -----------------------------------------------------

            # Starting pymol
        cmd.reinitialize()
        print("Loading structures to pymol...")
        structures = loading_structures_to_pymol(structure_files,query,cmd)
        cmd.remove("hydrogens")
        print("Superimposing structures...")
        structures = super_impose_structures(structures, max_initial_rmsd, cmd)
        print("Scoring structures based on similarity...")
        structures, align = SIMalign_algorithm(structures, max_dist)
        tmp_alignment_file_name = "User_data/"+job_key+"/results/alignment.aln"
        print("Updating sequence alignment...") # Is a bit slow. Could be optimized
        update_alignment_in_pymol(structures, align)
        cmd.save(tmp_alignment_file_name, selection="aln") 
            
        # ------------------------------- Hotspot finding ----------------------------------
            
        # ---------------- Hotspot functions -----------------
        
        def get_close_AA(structures, align):
                """
                Get close AA list and conserved regions
                :param structures: list of structures
                :param align: alignment
                :return: list of close AA
                """
                
                def get_core(structures):
                    """
                    Get conserved core region
                    :param structures: list of structures
                    :return: list of AA in conserved region
                    """
       
                    def median_40_percent(data):
                        """
                        Getting 40 percent median of the data
                        :param data: list of data
                        :return: 40 percent median of the data
                        """
                        n = len(data)
                        remove_count = int(0.2 * n)
                        # Directly pass the sliced list to statistics.median()
                        return statistics.median(sorted(data)[:n-remove_count])
                    
                    def findSurfaceResidues(selection="all", cutoff=2.5, doShow=0, quiet=1):
                        """
                        DESCRIPTION
                        
                        Finds those residues on the surface of a protein
                        that have at least 'cutoff' exposed A**2 surface area.
                        
                        USAGE
                        
                        findSurfaceResidues [ selection, [ cutoff, [ doShow ]]]
                        
                        ARGUMENTS
                        
                        selection = string: object or selection in which to find exposed
                        residues {default: all}
                        
                        cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}
                        
                        RETURNS
                        
                        (list: (chain, resv ) )
                        A Python list of residue numbers corresponding
                        to those residues w/more exposure than the cutoff.
                        
                        """
                        
                        def findSurfaceAtoms(selection="all", cutoff=2.5, quiet=1):
                            """
                            DESCRIPTION

                            Finds those atoms on the surface of a protein
                            that have at least 'cutoff' exposed A**2 surface area.
                            
                            USAGE

                            findSurfaceAtoms [ selection, [ cutoff ]]
                            
                            SEE ALSO
                            
                            findSurfaceResidues
                            """
                            cutoff, quiet = float(cutoff), int(quiet)
                            
                            tmpObj = cmd.get_unused_name("_tmp")
                            cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)
                            
                            cmd.set("dot_solvent", 1, tmpObj)
                            cmd.get_area(selection=tmpObj, load_b=1)
                            
                            # threshold on what one considers an "exposed" atom (in A**2):
                            cmd.remove(tmpObj + " and b < " + str(cutoff))
                            
                            selName = cmd.get_unused_name("exposed_atm_")
                            cmd.select(selName, "(" + selection + ") in " + tmpObj)
                            
                            cmd.delete(tmpObj)
                            
                            if not quiet:
                                print("Exposed atoms are selected in: " + selName)
                            return selName
                        
                        cutoff, doShow, quiet = float(cutoff), int(doShow), int(quiet)
                        
                        selName = findSurfaceAtoms(selection, cutoff, quiet)
                        
                        exposed = set()
                        cmd.iterate(selName, "exposed.add((chain,resv))", space=locals())
                        
                        selNameRes = cmd.get_unused_name("exposed_res_")
                        cmd.select(selNameRes, "byres " + selName)

                        if not quiet:
                            print("Exposed residues are selected in: " + selNameRes)
                            
                        if doShow:
                            cmd.show_as("spheres", "(" + selection + ") and polymer")
                            cmd.color("white", selection)
                            cmd.color("yellow", selNameRes)
                            cmd.color("red", selName)

                        return exposed
           
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
                    SurfaceResidues = set()
                    for residue in findSurfaceResidues(structures[0].first_chain):
                            SurfaceResidues.add(residue[1])

                    core = conserved_region - SurfaceResidues
                    return list(core)

                def update_closeAA_list(close_AA_list, tmp_close_AAs, count, residue, resi_list, index, align_seq_query, atomsCA_query, align_seq_template, atomsCA_template):
                    close_AAs = set()

                    # Adding glycines if any before because glycines don't have any side chain atoms
                    while len(resi_list[:resi_list.index(residue)]) > count:
                        close_AA_list[index].append(set())  # Add empty set for glycines
                        count += 1
                        
                    for ele in tmp_close_AAs:
                        close_AAs.add(resi_to_index(ele, align_seq_query, atomsCA_query))  # Update close_AAs with the computed residue index
                    close_AAs.discard(resi_to_index(residue, align_seq_template, atomsCA_template))
                    close_AA_list[index].append(close_AAs)  # Append the non-empty set
                    count += 1  # Increase count only once after appending

                    return close_AA_list, count
                
                core = get_core(structures)
                query_atoms = cmd.get_model(structures[0].side_chains).atom
                query_kd = cKDTree([atom.coord for atom in query_atoms])
                core_index = [resi_to_index(residue, align[0], structures[0].model.atom) for residue in core]
                close_AA_list = []
                for j, struc in enumerate(structures[1:]):
                    resi_list = [index_to_resi(index, align[j+1], struc.model.atom) for index in core_index]
                    resi_set = set(resi_list)

                    close_AA_list.append([])
                    atoms = cmd.get_model(struc.side_chains).atom
                    kd = cKDTree([atom.coord for atom in atoms])

                    residue = 0
                    tmp_close_AAs = set()
                    count = 0

                    for a in atoms:  # Remember, glycines have no side chains
                        if residue != int(a.resi):
                            if residue in resi_set:
                                close_AA_list, count = update_closeAA_list(close_AA_list, tmp_close_AAs, count, residue, resi_list, j, align[0], structures[0].model.atom, align[j+1], struc.model.atom)
                                tmp_close_AAs = set()
                            residue = int(a.resi)
                        if residue in resi_set:
                            # Finding atoms in query structure within 4Å of target template atom
                            tmp_set = set([int(query_atoms[x].resi) for x in query_kd.query_ball_point(a.coord, 4)])
                            tmp_close_AAs.update(tmp_set)

                    if residue in resi_set:
                        close_AA_list, count = update_closeAA_list(close_AA_list, tmp_close_AAs, count, residue, resi_list, j, align[0], structures[0].model.atom, align[j+1], struc.model.atom)

                    # Add empty set if sequence ends with glycine
                    while len(resi_list) > count:
                        close_AA_list[j].append(set())
                        count += 1

                return close_AA_list, core, core_index


        def findSurfaceAtoms(selection="all", cutoff=2.5, quiet=1):
                """
                DESCRIPTION

                Finds those atoms on the surface of a protein
                that have at least 'cutoff' exposed A**2 surface area.
                
                USAGE

                findSurfaceAtoms [ selection, [ cutoff ]]
                
                SEE ALSO
                
                findSurfaceResidues
                """
                cutoff, quiet = float(cutoff), int(quiet)
                
                tmpObj = cmd.get_unused_name("_tmp")
                cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)
                
                cmd.set("dot_solvent", 1, tmpObj)
                cmd.get_area(selection=tmpObj, load_b=1)
                
                # threshold on what one considers an "exposed" atom (in A**2):
                cmd.remove(tmpObj + " and b < " + str(cutoff))
                
                selName = cmd.get_unused_name("exposed_atm_")
                cmd.select(selName, "(" + selection + ") in " + tmpObj)
                
                cmd.delete(tmpObj)
                
                if not quiet:
                    print("Exposed atoms are selected in: " + selName)
                return selName

        def findSurfaceResidues(selection="all", cutoff=2.5, doShow=0, quiet=1):
                """
                DESCRIPTION
                
                Finds those residues on the surface of a protein
                that have at least 'cutoff' exposed A**2 surface area.
                
                USAGE
                
                findSurfaceResidues [ selection, [ cutoff, [ doShow ]]]
                
                ARGUMENTS
                
                selection = string: object or selection in which to find exposed
                residues {default: all}
                
                cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}
                
                RETURNS
                
                (list: (chain, resv ) )
                A Python list of residue numbers corresponding
                to those residues w/more exposure than the cutoff.
                
                """
                
                def findSurfaceAtoms(selection="all", cutoff=2.5, quiet=1):
                    """
                    DESCRIPTION

                    Finds those atoms on the surface of a protein
                    that have at least 'cutoff' exposed A**2 surface area.
                    
                    USAGE

                    findSurfaceAtoms [ selection, [ cutoff ]]
                    
                    SEE ALSO
                    
                    findSurfaceResidues
                    """
                    cutoff, quiet = float(cutoff), int(quiet)
                    
                    tmpObj = cmd.get_unused_name("_tmp")
                    cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)
                    
                    cmd.set("dot_solvent", 1, tmpObj)
                    cmd.get_area(selection=tmpObj, load_b=1)
                    
                    # threshold on what one considers an "exposed" atom (in A**2):
                    cmd.remove(tmpObj + " and b < " + str(cutoff))
                    
                    selName = cmd.get_unused_name("exposed_atm_")
                    cmd.select(selName, "(" + selection + ") in " + tmpObj)
                    
                    cmd.delete(tmpObj)
                    
                    if not quiet:
                        print("Exposed atoms are selected in: " + selName)
                    return selName
                
                cutoff, doShow, quiet = float(cutoff), int(doShow), int(quiet)
                
                selName = findSurfaceAtoms(selection, cutoff, quiet)
                
                exposed = set()
                cmd.iterate(selName, "exposed.add((chain,resv))", space=locals())
                
                selNameRes = cmd.get_unused_name("exposed_res_")
                cmd.select(selNameRes, "byres " + selName)

                if not quiet:
                    print("Exposed residues are selected in: " + selNameRes)
                    
                if doShow:
                    cmd.show_as("spheres", "(" + selection + ") and polymer")
                    cmd.color("white", selection)
                    cmd.color("yellow", selNameRes)
                    cmd.color("red", selName)

                return exposed
            
        cmd.extend("findSurfaceAtoms", findSurfaceAtoms)
        cmd.extend("findSurfaceResidues", findSurfaceResidues)

        def finding_hotspots(closeAA_list, align, structures, core, core_index, mode=1):
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
                                    for closeAA in closeAA_list[j][core_index.index(ref_index)]:
                                        if not not_bigger_AA(seq[closeAA],align[0][closeAA]):
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
                                                for close_to_close in closeAA_list[j][core_index.index(resi_off[0])]:
                                                    if not not_bigger_AA(seq[close_to_close],align[0][close_to_close]):
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
        
        def print_hotspot(hotspot_list, job_key, structures, mode=1):
                def not_bigger_AA(ref_AA,target_AA):
                    """
                    Check if target AA is not bigger than ref AA
                    :param target_AA: target AA
                    :param ref_AA: ref AA
                    :return: True if target AA is bigger than ref AA, False otherwise
                    """
                    if ref_AA == target_AA:
                        return True
                    elif target_AA == "G":
                        return True
                    elif target_AA == "A" and ref_AA not in {"G","P"}:
                        return True
                    elif target_AA == "V" and ref_AA == "I":
                        return True
                    elif target_AA == "F" and ref_AA == "Y":
                        return True
                    elif target_AA == "L" and ref_AA in {"F","Y"}:
                        return True
                    elif target_AA == "S" and ref_AA in {"C","T"}:
                        return True
                    else:
                        return False
                def bigger_AA(ref_AA,target_AA):
                    """
                    Check if target AA is bigger than ref AA
                    :param target_AA: target AA
                    :param ref_AA: ref AA
                    :return: True if target AA is bigger than ref AA, False otherwise
                    """
                    if ref_AA == target_AA:
                        return False
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
                doc_name="User_data/"+job_key+"/results/hotspots_mode_"+str(mode)+".html"
                doc=open(doc_name,"w")
                
                doc.write("------------------------------------------------------------------------------")
                if mode == 1:
                    doc.write("<h1>Printing possible single mutations in "+structures[0].name+"</h1>")
                elif mode == 2:
                    doc.write("<h1>Printing possible double mutations in "+structures[0].name+"</h1>")
                doc.write("------------------------------------------------------------------------------")
                doc.write("<h1>Legend</h1>")
                doc.write("<p style='color:rgb(255,9,0);'>Red : Mutations into smaller AAs</p>")
                doc.write("<p style='color:rgb(0,100,0);'>Green: Mutations into larger AAs</p>")
                doc.write("<p style='color:rgb(214,178,32);'>Yellow: Undefined</p>")
                doc.write("------------------------------------------------------------------------------")
                printed = set()
                for hotspot in hotspot_list:
                    for i, mutation_list in enumerate(hotspot[1]):
                        if not isinstance(mutation_list[0], str):
                            wt_list = hotspot[0]
                            non_comparable = hotspot[2][i]
                            
                            all_lower = False
                            all_higher = False
                            wt_print = ""
                            mu_print = ""
                            wt_list_sorted=sorted(wt_list)
                            mutation_list_sorted=sorted(mutation_list, key=lambda x: (x[0] if isinstance(x[0], int) else float('inf')))

                            for wt, mu in zip(wt_list_sorted, mutation_list_sorted):
                                if not_bigger_AA(wt[1], mu[1]):
                                    all_lower = True
                                elif bigger_AA(wt[1], mu[1]):
                                    all_higher = True
                                wt_print += wt[1]+str(wt[0])+" "
                                mu_print += mu[1]+" "
                            if not str(wt_print)+mu[1] in printed:
                                printed.add(str(wt_print)+mu[1])
                                newlist=[]
                                for j, x in enumerate(mutation_list_sorted): 
                                    if x[1] == mu[1]:
                                        newlist.append (structures[j+1].name)
                                if all_lower:
                                    text="<p style='color:rgb(255,9,0);'>"
                                    #doc.write("<p style='color:rgb(255,9,0);'>"+str(wt_print)+"-> "+str(mu_print)+" Non-comparable: "+", ".join([str(x) for x in non_comparable]) + "; Structures with given mutation:" + ", ".join([str(x) for x in newlist]) + "</p>")
                                elif all_higher:
                                    text="<p style='color:rgb(0,100,0);'>"
                                    #doc.write("<p style='color:rgb(0,255,4);'>"+str(wt_print)+" -> "+str(mu_print)+" Non-comparable"+str(len(non_comparable))+": "+", ".join([str(x) for x in non_comparable])+"; Structures with given mutation:" + ", ".join([str(x) for x in newlist]) + "</p>")
                                else:
                                    text="<p style='color:rgb(214,178,32);'>"
                                    #doc.write("<p style='color:rgb(214,178,32);'>"+str(wt_print)+" -> "+str(mu_print)+" Non-comparable"+str(len(non_comparable))+": "+", ".join([str(x) for x in non_comparable])+"; Structures with given mutation:" + ", ".join([str(x) for x in newlist]) + "</p>")
                                text+=str(wt_print)+" -> "+str(mu_print)
                                if len(non_comparable)>0:
                                    text+="; Non-comparable: "+", ".join([str(x) for x in non_comparable])
                                text+="; Structures with given mutation: " + ", ".join([str(x) for x in newlist]) + "</p>"
                                
                                doc.write(text)
                doc.write("------------------------------------------------------------------------------")
                doc.close()

            # ----------------------------------------------------
            
        print("Finding hotspots...")
        align = AlignIO.read(tmp_alignment_file_name,"clustal")
        # Getting closeAA list
        closeAA_list, core, core_index = get_close_AA(structures, align)
            
        hotspot_list_double = finding_hotspots(closeAA_list, align, structures, core, core_index, mode=2)
            
        hotspot_list_single = finding_hotspots(closeAA_list, align, structures, core, core_index, mode=1)
            
        hotspot_list = hotspot_list_double + hotspot_list_single
            
        print_hotspot(hotspot_list_double, job_key, structures, mode=2)
        print_hotspot(hotspot_list_single, job_key, structures, mode=1)
        """# Downloading results
        
        This might take a while
        """
            
            #@title Pymol visualization
            
            #structure_format = "cartoon" #@param ["spheres-sticks","cartoon","spheres","sticks"] {type:"string"}
            #@markdown   - `structure_format` specify how the structure should be showed in pymol.
            #show_in_pymol = "entire_chain_A" #@param ["only_not_conserved_core","only_core", "only_not_conserved","entire_chain_A","everything"] {type:"string"}
            #@markdown   - `show_in_pymol` specify what part of the structures that will be shown in pymol.
            #color_by_element = True #@param {type: "boolean"}
            #@markdown   - If `color_by_element` is ON then atom will be colored by element in pymol.
            
            # ---------------------------- Visualization functions ----------------------


        def color_by_score(struc):
                """
                Coloring by similarity in pymol
                :param struc: structure
                """
                def color_by_number(number):
                    """
                    DESCRIPTION
                    
                    Rainbow from red to white
                    number 0: red
                    number 1: white
                    
                    DEPENDENCIES
                    
                    import numpy as np
                    """
                    if number >= 1:
                        return [1,1,1]
                    return [0.8+(number/5),number,number]
                score = struc.score_list
                for i, atom in enumerate(struc.model.atom):
                    cmd.set_color(str(round(score[i],2)), color_by_number(score[i]))
                    cmd.color(str(round(score[i],2)), f"resi {atom.resi} and {struc.first_chain}")

        def structure_formating(structure_format,show_in_pymol):
                """
                Structure formating in pymol
                :param structure_format: structure format
                :param show_in_pymol: what part of the structures that will be shown in pymol
                """

                cmd.hide("everything","not chain A")
                if structure_format == "cartoon":
                    cmd.show_as("cartoon","not HETATM")
                elif structure_format == "spheres":
                    cmd.show_as("spheres","not HETATM")
                else:
                    cmd.show_as("sticks","not HETATM")
                    if structure_format == "spheres-sticks":
                        cmd.show("spheres","not HETATM")
                        cmd.set("sphere_transparency", "0.7")
                if show_in_pymol != "everything":
                    cmd.hide("everything","not chain A")
                    if show_in_pymol != "entire_chain_A" and show_in_pymol != "only_not_conserved":
                        cmd.hide("everything","not conserved")
                        if show_in_pymol != "only_not_conserved_core":
                            cmd.hide("everthing","conserved")
                    elif show_in_pymol == "only_not_conserved":
                        cmd.hide("everything","conserved")

            

            # --------------------------------------------------

        print("Coloring and selecting in pymol...")

            #Select conserved of query structure
        print("pre-coloring")
        print(structures[0].first_chain+select_by_list(list(structures[0].conserved_AA)))
        cmd.select("conserved", structures[0].first_chain+select_by_list(list(structures[0].conserved_AA)))
        print("post-coloring")
            #Select hotspot in query structure
            
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
                        # print(resi_list)
                        selection_string += f"{structures[i].first_chain} and ("
                        for resi in resi_list:
                            # print(resi)
                            # if not isinstance(resi, str):
                            selection_string += f"resi {resi[0]} or "
                        selection_string = selection_string[:-4]+")) or ("
                selection_string = selection_string[:-5]
                cmd.select(f"hotspot_{str(count)}", selection_string)
                count += 1
                # cmd.select("hotspots", structures[0].first_chain+select_by_list(list(structures[0].hotspot.keys())))
            
        structure_formating("cartoon","everything")
        #cmd.run("pymol-open-source/modules/pymol/util.py")  # Ensure utility functions are available    
        print("\tColoring by similarity:")
        for struc in structures:
                print(f"\t\tColoring {struc.name}")
                color_by_score(struc)
        #cmd.util.cnc("all")
        cmd.hide("cgo", "aln")
        cmd.set("seq_view_label_mode", "1")
        cmd.set("antialias", "4")
        cmd.set("ray_trace_mode", "1")
        zip_file = "SIMAlign_"+job_key#+".zip"
        pymol_file_name = "pymol_output"
        #alignment_file_name = "alignment.aln"
        alignment_file_name = "User_data/"+job_key+"/results/alignment.aln"
        #os.system(f"mv {tmp_alignment_file_name} {os.path.join(result_dir,alignment_file_name)}")
        cmd.save(os.path.join(result_dir,pymol_file_name+".pse"))
        #os.system(f"zip -r {zip_file} {result_dir}")
        shutil.make_archive("User_data/"+job_key+"/"+zip_file, "zip" ,"User_data/"+job_key+"/results", ".")
        print("Job is finished!")
        #try:
        #  files.download(zip_file)
        #  print(f'{zip_file} downloaded')
        #except Exception as e:
        #  print(f'\033[1;35mAn error occurred while downloading {zip_file}: {e}\033[0m')
        
    
    """# Instructions

    Code was made using ChatGPT as tool


    **Quick run:**
    1.   Press "Runtime" -> "Run all".
    2.   A buttom saying "Choose Files" will appear. Press it and choose the structure that you want to analyse.

    \\

    **Homology search method**

    Be aware that you can either run using foldseek or a user specified approach as homology search method.
    Pros and cons

    \\


    **License**

    Remember all packages used

    \\

    **Bugs:**

    - If you encounter any bugs, please report the issue to https://github.com/sokrypton/ColabFold/issues

    \\

    **Citing this work:**

    If you use our model please cite:

    Blaabjerg, L.M., Kassem, M.M., Good, L.L., Jonsson, N., Cagiada, M., Johansson, K.E., Boomsma, W., Stein, A. and Lindorff-Larsen, K., 2022. *Rapid protein stability prediction using deep learning representations*, bioRxiv. (https://doi.org/10.1101/2022.07.14.500157)
    **bold text**

    ```
    @article{blaabjerg2022rapid,
      title={Rapid protein stability prediction using deep learning representations},
      author={Blaabjerg, Lasse M and Kassem, Maher M and Good, Lydia L and Jonsson, Nicolas and Cagiada, Matteo and Johansson, Kristoffer E and Boomsma, Wouter and Stein, Amelie and Lindorff-Larsen, Kresten},
      journal={bioRxiv},
      year={2022},
      publisher={Cold Spring Harbor Laboratory}
    }
    ```
    """