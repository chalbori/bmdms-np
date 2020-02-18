from rdkit.Chem import AllChem as Chem
from rdkit.Chem import inchi
from tools import tool_chemical


class Chemical:
    __slots__ = [
        "_c_id",
        "_name",
        "_note",
        "_scaffold",
        "_orbitrap_qc",
        "_qtof_qc",
        "_ttof_qc",
        "_mol",
        "_smiles",
        "_inchi",
        "_inchikey",
        "_molecular_weight",
        "_molecular_formula",
        "_orbitrap_valid",
        "_qtof_valid",
    ]

    def __init__(self):
        self._c_id = None
        self._name = None
        self._note = None
        self._scaffold = None
        self._orbitrap_qc = None
        self._qtof_qc = None
        self._ttof_qc = None
        self._mol = None
        self._smiles = None
        self._inchi = None
        self._inchikey = None
        self._molecular_weight = None
        self._molecular_formula = None
        self._orbitrap_valid = None
        self._qtof_valid = None

    def __repr__(self):
        return (
            "<compound(c_id={c_id}, "
            "name={name}, "
            "note={note}, "
            "scaffold={scaffold},"
            "orbitrap_qc={orbitrap_qc},"
            "qtof_qc={qtof_qc},"
            "ttof_qc={ttof_qc}, "
            "mw={molecular_weight}, "
            "mf={molecular_formula}, "
            "smiles={smiles}, "
            "inchi={inchi}, "
            "inchikey={inchikey}, "
            "orbitrap_valid={orbitrap_valid},"
            "qtof_valid={qtof_valid}".format(
                c_id=self._c_id,
                name=self._name,
                note=self._note,
                scaffold=self._scaffold,
                orbitrap_qc=self._orbitrap_qc,
                qtof_qc=self._qtof_qc,
                ttof_qc=self._ttof_qc,
                molecular_weight=self._molecular_weight,
                molecular_formula=self._molecular_formula,
                smiles=self._smiles,
                inchi=self._inchi,
                inchikey=self._inchikey,
                orbitrap_valid=self._orbitrap_valid,
                qtof_valid=self._qtof_valid,
            )
        )

    def set_computable(self):
        mol = tool_chemical.read_string("mol", self._mol)
        # molecular_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
        # molecular_weight = Descriptors.ExactMolWt(mol)
        self._smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
        self._inchi = inchi.MolToInchi(mol)
        self._inchikey = inchi.MolToInchiKey(mol)
        self._molecular_formula = Chem.CalcMolFormula(mol)
        self._molecular_weight = Chem.CalcExactMolWt(mol)

    @property
    def c_id(self):
        return self._c_id

    @c_id.setter
    def c_id(self, c_id):
        self._c_id = c_id

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def note(self):
        return self._note

    @note.setter
    def note(self, note):
        self._note = note

    @property
    def scaffold(self):
        return self._scaffold

    @scaffold.setter
    def scaffold(self, scaffold):
        self._scaffold = scaffold

    @property
    def orbitrap_qc(self):
        return self._orbitrap_qc

    @orbitrap_qc.setter
    def orbitrap_qc(self, orbitrap_qc):
        self._orbitrap_qc = orbitrap_qc

    @property
    def qtof_qc(self):
        return self._qtof_qc

    @qtof_qc.setter
    def qtof_qc(self, qtof_qc):
        self._qtof_qc = qtof_qc

    @property
    def ttof_qc(self):
        return self._ttof_qc

    @ttof_qc.setter
    def ttof_qc(self, ttof_qc):
        self._ttof_qc = ttof_qc

    @property
    def mol(self):
        return self._mol

    @mol.setter
    def mol(self, mol):
        self._mol = mol

    @property
    def smiles(self):
        return self._smiles

    @smiles.setter
    def smiles(self, smiles):
        self._smiles = smiles

    @property
    def inchi(self):
        return self._inchi

    @inchi.setter
    def inchi(self, inchi):
        self._inchi = inchi

    @property
    def inchikey(self):
        return self._inchikey

    @inchikey.setter
    def inchikey(self, inchikey):
        self._inchikey = inchikey

    @property
    def molecular_weight(self):
        return self._molecular_weight

    @molecular_weight.setter
    def molecular_weight(self, molecular_weight):
        self._molecular_weight = molecular_weight

    @property
    def molecular_formula(self):
        return self._molecular_formula

    @molecular_formula.setter
    def molecular_formula(self, molecular_formula):
        self._molecular_formula = molecular_formula

    @property
    def orbitrap_valid(self):
        return self._orbitrap_valid

    @orbitrap_valid.setter
    def orbitrap_valid(self, orbitrap_valid):
        self._orbitrap_valid = orbitrap_valid

    @property
    def qtof_valid(self):
        return self._qtof_valid

    @qtof_valid.setter
    def qtof_valid(self, qtof_valid):
        self._qtof_valid = qtof_valid
