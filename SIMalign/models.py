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
        """Get fasta sequence from the structure."""
        lines = self.cmd.get_fastastr(self.name).split("\n")
        fasta = ""
        for line in lines[1:]:
            if line.startswith(">"):
                return fasta
            else:
                fasta += line
        return fasta

    def validate_structure_format(self):
        """Validate the structure format and check for non-canonical amino acids."""
        self.stored.residues = []
        self.cmd.iterate(self.name, "stored.residues.append(resi)")
        if [resi for resi in self.stored.residues if not resi.isdigit()]:
            error_msg = f"{self.name} was removed due to invalid residues."
            print(f"\033[35m\t{error_msg}\033[0m")
            self.cmd.delete(self.name)
            return error_msg

        if len(self.get_fasta()) != len(self.cmd.get_model(self.CA).atom):
            error_msg = f"{self.name} was removed because it contains non-canonical amino acids."
            print(f"\033[35m\t{error_msg}\033[0m")
            self.cmd.delete(self.name)
            return error_msg

        return None  # Validation passed
