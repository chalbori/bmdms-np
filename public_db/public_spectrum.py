from decimal import Decimal
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import inchi


def _num_decimal_places(formatted_peaklist):
    # type check. formatted_peaklist must be a list of list(str)
    if not issubclass(type(formatted_peaklist[0][0]), str):
        raise TypeError(
            "{!r} is not a subclass of {!r}".format(type(formatted_peaklist), str)
        )

    min_dp = 5
    dp_cnt_generator = (
        Decimal(peak[0]).as_tuple().exponent for peak in formatted_peaklist
    )
    dp_cnt_dict = dict()
    all_dp_cnt = 0
    for dp in dp_cnt_generator:
        all_dp_cnt += 1
        each_dp_cnt = dp_cnt_dict.get(dp, None)
        if each_dp_cnt:
            dp_cnt_dict[each_dp_cnt] += 1
        else:
            dp_cnt_dict.update({dp: 0})

    dp_percent_dict = {dp: dp_cnt_dict.get(dp) / all_dp_cnt for dp in dp_cnt_dict}

    for dp in dp_percent_dict.keys():
        if dp_percent_dict.get(dp) >= 0.8:
            return dp
        elif dp < min_dp:
            min_dp = dp
    return min_dp


def _convert_into_mol(structure_info, info_type="mol"):
    if info_type == "mol block":
        return Chem.MolFromMolBlock(structure_info)
    elif info_type == "smiles":
        return Chem.MolFromSmiles(structure_info)
    else:
        return structure_info


class SpectrumError(Exception):
    pass


class Spectrum:
    __slots__ = [
        "_spectrum_id",
        "_precursor_type",
        "_precursor_mz",
        "_instrument_type",
        "_instrument_name",
        "_collision_energy",
        "_author",
        "_peaklist",
        "_compound_name",
        "_molecular_weight",
        "_molecular_formula",
        "_inchi",
        "_inchikey",
        "_smiles",
        "_tags",
        "_num_decimal_places",
    ]

    @classmethod
    def from_exp_record(
        cls,
        precursor_type=None,
        precursor_mz=None,
        instrument_type=None,
        instrument_name=None,
        collision_energy=None,
        author=None,
        spectrum_id=None,
        peaklist=None,
        compound_name=None,
        tags=None,
        structure_str=None,
        str_type=None,
    ):

        spec = cls()

        try:  # set num_decimal_places.
            spec.num_decimal_places = _num_decimal_places(peaklist)
        except ValueError as e:
            raise SpectrumError(e.args) from e

        try:
            spec.precursor_type = precursor_type
            spec.precursor_mz = precursor_mz
            spec.instrument_name = instrument_name
            spec.instrument_type = instrument_type
            spec.collision_energy = collision_energy
            spec.author = author
            spec.spectrum_id = spectrum_id
            spec.peaklist = [[float(x[0]), float(x[1])] for x in peaklist]
            spec.compound_name = compound_name
            spec.tags = tags
        except (TypeError, ValueError) as e:
            raise SpectrumError(
                "Spectrum Initiation Error:\n{!r}".format(e.args)
            ) from e
        except Exception as e:
            print(e)

        # set computables
        try:
            spec.set_computables(structure_str=structure_str, str_type=str_type)

        except SpectrumError as e:
            print(e)

        # if spec.check_integrity():
        #     raise AssertionError()
        return spec

    def __init__(self):
        self._spectrum_id = None
        self._precursor_type = None
        self._precursor_mz = None
        self._instrument_type = None
        self._instrument_name = None
        self._collision_energy = None
        self._author = None
        self._peaklist = None
        self._compound_name = None
        self._molecular_weight = None
        self._molecular_formula = None
        self._inchi = None  # TEXT
        self._inchikey = None  # TEXT
        self._smiles = None  # TEXT
        self._tags = None  # TEXT
        self._num_decimal_places = None

    def set_computables_from_mol(self, mol):
        try:  # warning comes up in pycharm (bug of pycharm)
            self.molecular_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
            self.molecular_weight = Descriptors.ExactMolWt(mol)
            self.inchi = inchi.MolToInchi(mol)
            self.inchikey = inchi.MolToInchiKey(mol)
            self.smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
        except Exception as e:
            raise SpectrumError(
                "Error occurred while computing properties" + e.args
            ) from e

        assert self.molecular_formula is not None, "molecular-formula can't be None"
        assert self.molecular_weight is not None, "molecular-weight can't be None"
        assert self.inchi is not None, "inchi can't be None"
        assert self.inchikey is not None, "inchikey can't be None"
        assert self.smiles is not None, "smiles can't be None"

    def set_computables(self, structure_str, str_type="mol"):
        if not issubclass(type(structure_str), str):
            raise TypeError(
                "mol generation error: {!r} is not supported.\n Supported type: {!r}".format(
                    type(structure_str), str
                )
            )
        if not structure_str:  # True if string is empty
            raise ValueError(
                "mol generation error : {!r} is an empty string.".format(structure_str)
            )
        if not structure_str.strip():  # True if string is whitespace-only
            raise ValueError(
                "mol generation error: {!r} is an whitespace-only string".format(
                    structure_str
                )
            )
        _mol = _convert_into_mol(structure_str, info_type=str_type)
        self.set_computables_from_mol(_mol)

        # self.mass_error = _calc_mass_error(self.molecular_formula, self.precursor_type, self.precursor_mz)

    @property
    def spectrum_id(self):
        return self._spectrum_id

    @spectrum_id.setter
    def spectrum_id(self, spectrum_id):
        self._spectrum_id = spectrum_id

    @property
    def precursor_type(self):
        return self._precursor_type

    @precursor_type.setter
    def precursor_type(self, precursor_type):
        self._precursor_type = precursor_type

    @property
    def precursor_mz(self):
        return self._precursor_mz

    @precursor_mz.setter
    def precursor_mz(self, precursor_mz):
        self._precursor_mz = precursor_mz

    @property
    def instrument_type(self):
        return self._instrument_type

    @instrument_type.setter
    def instrument_type(self, instrument_type):
        self._instrument_type = instrument_type

    @property
    def instrument_name(self):
        return self._instrument_name

    @instrument_name.setter
    def instrument_name(self, instrument_name):
        self._instrument_name = instrument_name

    @property
    def collision_energy(self):
        return self._collision_energy

    @collision_energy.setter
    def collision_energy(self, collision_energy):
        self._collision_energy = collision_energy

    @property
    def author(self):
        return self._author

    @author.setter
    def author(self, author):
        self._author = author

    @property
    def peaklist(self):
        return self._peaklist

    @peaklist.setter
    def peaklist(self, peaklist):
        self._peaklist = peaklist

    @property
    def compound_name(self):
        return self._compound_name

    @compound_name.setter
    def compound_name(self, compound_name):
        self._compound_name = compound_name

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, tags):
        self._tags = tags

    @property
    def molecular_formula(self):
        return self._molecular_formula

    @molecular_formula.setter
    def molecular_formula(self, molecular_formula):
        self._molecular_formula = molecular_formula

    @property
    def molecular_weight(self):
        return self._molecular_weight

    @molecular_weight.setter
    def molecular_weight(self, molecular_weight):
        self._molecular_weight = molecular_weight

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
    def smiles(self):
        return self._smiles

    @smiles.setter
    def smiles(self, smiles):
        self._smiles = smiles

    @property
    def num_decimal_places(self):
        return self._num_decimal_places

    @num_decimal_places.setter
    def num_decimal_places(self, num_decimal_places):
        self._num_decimal_places = num_decimal_places
