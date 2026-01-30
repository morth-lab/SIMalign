from scipy.spatial import cKDTree


def threeletter2oneletter(resn: str) -> str:

    AA3_TO_AA1 = {
        # Canonical 20
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
        "GLU":"E","GLN":"Q","GLY":"G","HIS":"H","ILE":"I",
        "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
        "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",

        # Common non-canonical / genetically encoded
        "SEC":"U",  # selenocysteine
        "PYL":"O",  # pyrrolysine

        # Very common crystallography substitution
        "MSE":"M",  # selenomethionine -> Met

        # Histidine naming variants (Amber/CHARMM/protonation states)
        "HID":"H","HIE":"H","HIP":"H",
        "HSD":"H","HSE":"H","HSP":"H",

        # Protonated acids (common in some forcefields)
        "ASH":"D",  # protonated Asp
        "GLH":"E",  # protonated Glu

        # Cys variants (forcefield / disulfide annotations)
        "CYX":"C",  # disulfide-bonded cysteine often labeled CYX
        "CYM":"C",  # deprotonated cysteine

        # Lys variants
        "LYN":"K",  # neutral lysine (deprotonated)

        # Common phosphorylated residues (often used as residue names)
        "SEP":"S",  # phosphoserine
        "TPO":"T",  # phosphothreonine
        "PTR":"Y",  # phosphotyrosine

        # A few other common modified AA-like residues you may meet
        "HYP":"P",  # hydroxyproline
        "CSO":"C",  # cysteine sulfinic acid (often treated as C)
        "CSD":"C",  # cysteine sulfinic/sulfonic variants (often treated as C)
        "CSS":"C",  # cysteine-related variant
    }

    out = AA3_TO_AA1.get(resn.upper(), "X")  # X for unknown
    if out == "X":
        print(f'<p style="color:orange;"><b>Warning:</b> Unknown residue name encountered: {resn}. Mapped to "X".</p>')
    return out



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
        if cmd.get_model(self.not_HETATM).atom[0].chain:
            self.first_chain = f"{self.not_HETATM} and chain {cmd.get_model(self.not_HETATM).atom[0].chain}"
        else:

            print(f'<p style="color:orange;"><b>Warning:</b> No chain information found for {self.name}, so the whole structure will be used.</p>')
            self.first_chain = self.not_HETATM

        self.CA = f"{self.first_chain} and name CA"
        self.side_chains = f"{self.first_chain} and (not name CA and not name N and not name C and not name O)"
        self.model = None
        self.cKDTree = None
        self.score_list = []
        self.conserved_AA = set()
    
    def get_fasta(self):
        fasta = ""
        for atom in self.model.atom:
            fasta += threeletter2oneletter(atom.resn)
        return fasta

    def pymol_build_in_get_fasta(self):
        """Get fasta sequence from the structure."""
        lines = self.cmd.get_fastastr(self.first_chain).split("\n")
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
            print(f'<p style="color:red;"><b>ERROR:</b> {error_msg}</p>')
            self.cmd.delete(self.name)
            return error_msg
        


        if len(self.pymol_build_in_get_fasta()) != len(self.cmd.get_model(self.CA).atom):
            error_msg = f"{self.name} was removed because the FASTA and CA atom counts do not match."
            print(f'<p style="color:red;"><b>ERROR:</b> {error_msg}</p>')
            fasta_length = len(self.pymol_build_in_get_fasta())
            ca_length = len(self.cmd.get_model(self.CA).atom)
            print(f'<p style="color:red;"><b>DETAILS:</b> Fasta length: {fasta_length}, CA atom length: {ca_length}</p>')
            print(self.cmd.get_model(self.not_HETATM).atom[0].chain)
            print("build in FASTA:")
            print(self.pymol_build_in_get_fasta())
            print("custom FASTA:")
            print(self.get_fasta())
            print("ZZZ")
            self.cmd.delete(self.name)
            return error_msg

        return None  # Validation passed
