from scipy.spatial import cKDTree

class StructureFile:
    """Represents a structure file with its name, location, and type."""
    def __init__(self, name, file_location, type):
        self.name = name
        self.file_location = file_location
        self.type = type


class Structure:
    """Represents a structure loaded into PyMOL with associated attributes."""
    def __init__(self, name, cmd, stored):
        self.name = name
        self.cmd = cmd
        self.stored = stored
        self.not_HETATM = f"{self.name} and not HETATM"
        self.first_chain = f"{self.not_HETATM} and chain {cmd.get_model(self.not_HETATM).atom[0].chain}"
        self.CA = f"{self.first_chain} and name CA"
        self.side_chains = f"{self.first_chain} and (not name CA and not name N and not name C and not name O)"
        self.model = None
        self.cKDTree = None
        self.score_list = []
        self.conserved_AA = set()
    
    def get_fasta(self):
        """
        Get fasta sequence from structure
        """
        lines = self.cmd.get_fastastr(self.name).split("\n")
        fasta = ""
        for line in lines[1:]:
            if line.startswith(">"):
                return fasta
            else:
                fasta += line
        return fasta

    def test_struc_format(self):
        """Check if the structure format is valid."""
        self.stored.residues = []
        self.cmd.iterate(self.name,"stored.residues.append(resi)")
        if [resi for resi in self.stored.residues if not resi.isdigit()]:
            print("\033[35m\t"+self.name,"was removed due to error in its residues\033[0m")
            self.cmd.delete(self.name)
        if len(self.get_fasta()) != len(self.cmd.get_model(self.CA).atom):
            print("\033[35m\t"+self.name,"was removed because it contain non canonical amino acids\033[0m")
            self.cmd.delete(self.name)

    # def test_struc_format(self):
    #     """Check if the structure format is valid."""
    #     stored_residues = []
    #     self.cmd.iterate(self.name, "stored.residues.append(resi)")
    #     if any(not resi.isdigit() for resi in stored_residues):
    #         print(f"\033[35m{self.name} removed due to invalid residues.\033[0m")
    #         self.cmd.delete(self.name)

    # def calculate_score(self, ref_structure, max_dist):
    #     """Calculate similarity scores for a structure compared to a reference structure."""
    #     self.model = self.cmd.get_model(self.CA)
    #     self.cKDTree = cKDTree([atom.coord for atom in self.model.atom])
    #     ref_kd = ref_structure.cKDTree
    #     scores = [
    #         sum(1.0 / (1 + ref_kd.query(atom.coord)[0]) for atom in self.model.atom)
    #         for atom in ref_structure.model.atom
    #     ]
    #     self.score_list = scores