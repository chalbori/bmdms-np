import sqlalchemy
from sqlalchemy import Column, Integer, Float
from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy import or_
from rdkit.Chem.rdmolfiles import SDWriter
from data.chemical import Chemical
from tools import tool_chemical as tc

Base = declarative_base()


class CompoundInfo(Base):
    __tablename__ = "info_compound"
    c_id = Column(mysql.VARCHAR(12), primary_key=True)
    name = Column(mysql.VARCHAR(1000))
    note = Column(mysql.VARCHAR(30))
    scaffold = Column(mysql.VARCHAR(30))
    orbitrap_qc = Column(mysql.VARCHAR(10))
    qtof_qc = Column(mysql.VARCHAR(10))
    ttof_qc = Column(mysql.VARCHAR(10))
    mol = Column(mysql.LONGTEXT)
    smiles = Column(mysql.TEXT)
    inchi = Column(mysql.TEXT)
    inchikey = Column(mysql.VARCHAR(30))
    molecular_weight = Column(Float)
    molecular_formula = Column(mysql.VARCHAR(200))
    orbitrap_valid = Column(Integer)
    qtof_valid = Column(Integer)

    @staticmethod
    def from_chemical(chemical):
        _chem_info = CompoundInfo()
        _chem_info.c_id = chemical.c_id
        _chem_info.name = chemical.name
        _chem_info.note = chemical.note
        _chem_info.scaffold = chemical.scaffold
        _chem_info.orbitrap_qc = chemical.orbitrap_qc
        _chem_info.qtof_qc = chemical.qtof_qc
        _chem_info.ttof_qc = chemical.ttof_qc
        _chem_info.mol = chemical.mol
        _chem_info.smiles = chemical.smiles
        _chem_info.inchi = chemical.inchi
        _chem_info.inchikey = chemical.inchikey
        _chem_info.molecular_weight = chemical.molecular_weight
        _chem_info.molecular_formula = chemical.molecular_formula
        _chem_info.orbitrap_valid = chemical.orbitrap_valid
        _chem_info.qtof_valid = chemical.qtof_valid
        return _chem_info

    def as_chemical(self):
        _chemical = Chemical()
        _chemical.c_id = self.c_id
        _chemical.name = self.name
        _chemical.note = self.note
        _chemical.scaffold = self.scaffold
        _chemical.orbitrap_qc = self.orbitrap_qc
        _chemical.qtof_qc = self.qtof_qc
        _chemical.ttof_qc = self.ttof_qc
        _chemical.mol = self.mol
        _chemical.smiles = self.smiles
        _chemical.inchi = self.inchi
        _chemical.inchikey = self.inchikey
        _chemical.molecular_weight = self.molecular_weight
        _chemical.molecular_formula = self.molecular_formula
        _chemical.orbitrap_valid = self.orbitrap_valid
        _chemical.qtof_valid = self.qtof_valid
        return _chemical

    def __repr__(self):
        return (
            "<compound_info(c_id={c_id}, "
            "name={name}, "
            "note={note}, "
            "scaffold={scaffold}, "
            "mw={molecular_weight}, "
            "mf={molecular_formula}, "
            "smiles={smiles}, "
            "inchi={inchi}, "
            "inchikey={inchikey}, "
            "orbitrap_valid={orbitrap_valid}, "
            "qtof_valid={qtof_valid}".format(
                c_id=self.c_id,
                name=self.name,
                note=self.note,
                scaffold=self.scaffold,
                molecular_weight=self.molecular_weight,
                molecular_formula=self.molecular_formula,
                orbitrap_qc=self.orbitrap_qc,
                qtof_qc=self.qtof_qc,
                ttof_qc=self.ttof_qc,
                smiles=self.smiles,
                inchi=self.inchi,
                inchikey=self.inchikey,
                orbitrap_valid=self.orbitrap_valid,
                qtof_valid=self.qtof_valid,
            )
        )


class CompoundRepository:
    def __init__(self, db_uri, echo=False):
        engine = sqlalchemy.create_engine(db_uri, echo=echo)
        Base.metadata.create_all(engine)
        Session = sessionmaker()
        Session.configure(bind=engine)
        self.session = Session()

    def get_stream(self):
        query = self.session.query(CompoundInfo).order_by(CompoundInfo.c_id)
        for result in query.all():
            yield result.as_chemical()
        
    def get_stream_valid(self):
        query = (
            self.session.query(CompoundInfo)
            .filter(or_(CompoundInfo.orbitrap_valid == 1, CompoundInfo.qtof_valid == 1))
            .order_by(CompoundInfo.c_id)
        )
        for result in query.all():
            yield result

    def get_orbitrap_valid_stream(self, orbitrap_valid=1):
        query = (
            self.session.query(CompoundInfo)
            .filter(CompoundInfo.orbitrap_valid == orbitrap_valid)
            .order_by(CompoundInfo.c_id)
        )
        for result in query.all():
            yield result

    def get_qtof_valid_stream(self, qtof_valid=1):
        query = (
            self.session.query(CompoundInfo)
            .filter(CompoundInfo.qtof_valid == qtof_valid)
            .order_by(CompoundInfo.c_id)
        )
        for result in query.all():
            yield result

    def get(self, c_id="%"):
        query = self.session.query(CompoundInfo).filter(CompoundInfo.c_id.like(c_id))
        return query.first().as_chemical()

    def get_list_by(self, c_id="%"):
        result = []
        query = self.session.query(CompoundInfo).filter(CompoundInfo.c_id.like(c_id))
        for record in query.all():
            result.append(record.as_chemical())
        return result

    def get_list_inchikey(self):
        result = []
        query = (
            self.session.query(CompoundInfo.inchikey)
            .distinct()
            .filter(or_(CompoundInfo.orbitrap_valid == 1, CompoundInfo.qtof_valid == 1))
            .order_by(CompoundInfo.inchikey)
        )
        for record in query.all():
            result.append(record[0])
        return result

    def get_list_smiles_metfrag(self):
        result = []
        query = (
            self.session.query(CompoundInfo.smiles)
            .distinct()
            .filter(or_(CompoundInfo.orbitrap_valid == 1, CompoundInfo.qtof_valid == 1))
            .order_by(CompoundInfo.smiles)
        )
        for record in query.all():
            result.append(record[0])
        return result

    def get_c_id_list(self, orbitrap_valid=1, qtof_valid=1):
        result = []
        query = self.session.query(CompoundInfo.c_id).filter(
            or_(
                CompoundInfo.orbitrap_valid == orbitrap_valid,
                CompoundInfo.qtof_valid == qtof_valid,
            )
        )
        for record in query.all():
            result.append(record.c_id)
        return result

    def update_identifier(self, dict_identifier):
        """[summary]
        
        Arguments:
            dict_identifier {[type]}
             -- {compound_id:{'smiles':smiles,'inchi':smiles,'inchikey':inchikey}}
        """
        for key in dict_identifier.keys():
            self.session.query(CompoundInfo).filter(CompoundInfo.c_id == key).update(
                {
                    CompoundInfo.smiles: dict_identifier[key]["smiles"],
                    CompoundInfo.inchi: dict_identifier[key]["inchi"],
                    CompoundInfo.inchikey: dict_identifier[key]["inchikey"],
                    CompoundInfo.molecular_weight: dict_identifier[key][
                        "molecular_weight"
                    ],
                    CompoundInfo.molecular_formula: dict_identifier[key][
                        "molecular_formula"
                    ],
                }
            )
            self.session.flush()
        self.session.commit()

    def write_sdf_file(self, c_id="%", output_file=""):
        list_compounds = self.get_list_by(c_id=c_id)
        writer = SDWriter(output_file)
        for compound in list_compounds:
            try:
                mol = tc.read_string("mol", compound.mol)
                mol.SetProp("_Name", compound.name)
                mol.SetProp("_Compound_id", compound.c_id)
                mol.SetProp("_Scaffold", compound.scaffold if compound.scaffold else "")
                mol.SetProp("_SMILES", compound.smiles)
                mol.SetProp("_INCHI", compound.inchi)
                mol.SetProp("_INCHIKEY", compound.inchikey)
                mol.SetProp("_Molecular_weight", str(compound.molecular_weight))
                mol.SetProp("_Molecular_formula", compound.molecular_formula)
                writer.write(mol)
            except:
                print("c_id: {} has an error!".format(compound.c_id))

    def add(self, chemical_list):
        for chemical in chemical_list:
            if not issubclass(type(chemical), Chemical):
                raise TypeError(
                    "{!r} is an unsupported type\n"
                    "{!r} is expected.".format(type(chemical), Chemical)
                )
            else:
                pass

        wrapped_chemical_list = []
        for chemical in chemical_list:
            wrapped_chemical = CompoundInfo.from_chemical(chemical=chemical)
            wrapped_chemical_list.append(wrapped_chemical)
        self.session.bulk_save_objects(wrapped_chemical_list)
        self.session.commit()

    def close(self):
        self.session.close()