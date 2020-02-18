
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import inchi


def convert(input, input_mod='smi'):
    """
    convert SMILES into other molecular identifier
    :param input: SMILES
    :param input_mod: 'smi'
    :return: str(molecular formula), str(inchi), str(inchikey)
    """
    mol = read_string(input_mod, input)
    molecular_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    molecular_inchi = inchi.MolToInchi(mol)
    molecular_inchikey = inchi.MolToInchiKey(mol)
    return molecular_formula, molecular_inchi, molecular_inchikey


def calculate_mw(input, input_mod='smi'):
    """
    calculate exact molecular weight
    :param input:
    :param input_mod:
    :return:
    """
    mol = read_string(input_mod, input)
    molecular_weight = Descriptors.ExactMolWt(mol)
    return molecular_weight

def read_string(mode, string, **kwargs):
    """Read in a molecule from a string.
    Required parameters:
       mode - see the informats variable for a list of available
                input formats
       string
    Example:
    >>> input = "C1=CC=CS1"
    >>> mymol = read_string("smi", input)
    >>> len(mymol.atoms)
    """
    mode = mode.lower()
    if mode == "mol" or mode == "sdf":
        mol = Chem.MolFromMolBlock(string, **kwargs)
    elif mode == "mol2":
        mol = Chem.MolFromMol2Block(string, **kwargs)
    elif mode == "pdb":
        mol = Chem.MolFromPDBBlock(string, **kwargs)
    elif mode == "smi":
        s = string.split()
        mol = Chem.MolFromSmiles(s[0], **kwargs)
        if mol:
            mol.SetProp("_Name", ' '.join(s[1:]))
    elif mode == 'inchi' and Chem.INCHI_AVAILABLE:
        mol = Chem.inchi.MolFromInchi(string, **kwargs)
    else:
        raise ValueError("%s is not a recognised RDKit mode" % mode)
    # return Molecule(mol)
    return mol
