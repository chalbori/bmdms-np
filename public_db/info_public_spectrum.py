import ast
import sqlalchemy
from sqlalchemy import Column, String, Integer, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import func
from public_db.public_spectrum import Spectrum

BASE = declarative_base()


class SpecInfo(BASE):
    __tablename__ = "info"

    spectrum_id = Column(String, primary_key=True)
    precursor_type = Column(String)
    precursor_mz = Column(Float)
    instrument_type = Column(String)
    instrument_name = Column(String)
    collision_energy = Column(String)
    author = Column(String)
    peaklist = Column(String)
    compound_name = Column(String)
    molecular_weight = Column(Float)
    molecular_formula = Column(String)
    inchi = Column(String)
    inchikey = Column(String)
    smiles = Column(String)
    tags = Column(String)
    num_decimal_places = Column(Integer)

    def __repr__(self):
        return (
            "<Info(spectrum_id={spectrum_id}, "
            "precursor_type={precursor_type}, "
            "precursor_mz={precursor_mz}, "
            "instrument_type={instrument_type}, "
            "instrument_name={instrument_name}, "
            "collision_energy={collision_energy}, "
            "author={author}, "
            "peaklist={peaklist}, "
            "compound_name={compound_name}, "
            "molecular_weight={molecular_weight}, "
            "molecular_formula={molecular_formula}, "
            "inchi={inchi}, "
            "inchikey={inchikey}, "
            "smiles={smiles}, "
            "tags={tags}, "
            "num_decimal_places={num_decimal_places}".format(
                spectrum_id=self.spectrum_id,
                precursor_type=self.precursor_type,
                precursor_mz=self.precursor_mz,
                instrument_type=self.instrument_type,
                instrument_name=self.instrument_name,
                collision_energy=self.collision_energy,
                author=self.author,
                peaklist=self.peaklist,
                compound_name=self.compound_name,
                molecular_weight=self.molecular_weight,
                molecular_formula=self.molecular_formula,
                inchi=self.inchi,
                inchikey=self.inchikey,
                smiles=self.smiles,
                tags=self.tags,
                num_decimal_places=self.num_decimal_places,
            )
        )

    def as_spectrum(self):
        _spectrum = Spectrum()
        _spectrum.spectrum_id = self.spectrum_id
        _spectrum.precursor_type = self.precursor_type
        _spectrum.precursor_mz = self.precursor_mz
        _spectrum.instrument_type = self.instrument_type
        _spectrum.instrument_name = self.instrument_name
        _spectrum.collision_energy = self.collision_energy
        _spectrum.author = self.author
        _spectrum.peaklist = [
            [float(x[0]), float(x[1])] for x in ast.literal_eval(self.peaklist)
        ]
        _spectrum.compound_name = self.compound_name
        _spectrum.molecular_weight = self.molecular_weight
        _spectrum.molecular_formula = self.molecular_formula
        _spectrum.inchi = self.inchi
        _spectrum.inchikey = self.inchikey
        _spectrum.smiles = self.smiles
        _spectrum.tags = self.tags
        _spectrum.num_decimal_places = self.num_decimal_places
        return _spectrum


class LocalSpectrumRepository:
    def __init__(self, db_uri):
        engine = sqlalchemy.create_engine(db_uri, echo=False)
        BASE.metadata.create_all(engine)
        Session = sessionmaker()
        Session.configure(bind=engine)
        self.session = Session()

    def get_smiles_stream(self, instrument_type="%", precursor_type="%", tags="%"):
        query = (
            self.session.query(SpecInfo.smiles)
            .distinct(SpecInfo.smiles)
            .filter(
                SpecInfo.instrument_type.like(instrument_type),
                SpecInfo.precursor_type.like(precursor_type),
                SpecInfo.tags.like(tags),
            )
        )
        for result in query.all():
            yield result.smiles  # smiles str.

    def get_smiles_list(self, instrument_type="%", precursor_type="%", tags="%"):
        list_return = []
        query = (
            self.session.query(SpecInfo.smiles)
            .distinct(SpecInfo.smiles)
            .filter(
                SpecInfo.instrument_type.like(instrument_type),
                SpecInfo.precursor_type.like(precursor_type),
                SpecInfo.tags.like(tags),
            )
        )
        for result in query.all():
            list_return.append(result.smiles)  # smiles str.
        return list_return

    def get_inchikey_list(self, instrument_type="%", precursor_type="%", tags="%"):
        list_return = []
        query = (
            self.session.query(SpecInfo.inchikey)
            .distinct(SpecInfo.inchikey)
            .filter(
                SpecInfo.instrument_type.like(instrument_type),
                SpecInfo.precursor_type.like(precursor_type),
                SpecInfo.tags.like(tags),
            )
        )
        for result in query.all():
            list_return.append(result.inchikey)  # inchikey str.
        return list_return

    def get_spectrum_dict_by_inchikey_choose(
        self,
        inchikey_list,
        instrument_type="%",
        precursor_type="%",
        include_len_0=True
    ):
        spectrum_inchikey_dict = dict()
        total_size = len(spectrum_inchikey_dict.keys())
        for n, inchikey in enumerate(inchikey_list):
            # print('Making target spectrum dict : {} / {}'.format(n+1, len(inchikey_list)))
            # print("get spectrum dict one smiles!")
            spectrum_list = []
            query = self.session.query(SpecInfo).filter(
                SpecInfo.inchikey == inchikey,
                SpecInfo.instrument_type == instrument_type,
                SpecInfo.precursor_type == precursor_type
            )
            for record in query.all():
                result_spectrum = record.as_spectrum()
                if len(result_spectrum.peaklist) > 4:
                    spectrum_list.append(record.as_spectrum())
            if include_len_0:
                spectrum_inchikey_dict[inchikey] = spectrum_list
                # spectrum_smiles_dict[smiles].append(spectrum_list)
            elif include_len_0 is False:
                if len(spectrum_list):
                    spectrum_inchikey_dict[inchikey] = spectrum_list
                    # spectrum_smiles_dict[smiles].append(spectrum_list)
        return spectrum_inchikey_dict

    def get_spectrum_dict_by_inchikey_choose_tags(
        self,
        inchikey_list,
        instrument_type="%",
        precursor_type="%",
        tags="%",
        include_len_0=True
    ):
        spectrum_inchikey_dict = dict()
        total_size = len(spectrum_inchikey_dict.keys())
        for n, inchikey in enumerate(inchikey_list):
            # print('Making target spectrum dict : {} / {}'.format(n+1, len(inchikey_list)))
            # print("get spectrum dict one smiles!")
            spectrum_list = []
            query = self.session.query(SpecInfo).filter(
                SpecInfo.inchikey == inchikey,
                SpecInfo.instrument_type == instrument_type,
                SpecInfo.precursor_type == precursor_type,
                SpecInfo.tags == tags
            )
            for record in query.all():
                result_spectrum = record.as_spectrum()
                if len(result_spectrum.peaklist) > 4:
                    spectrum_list.append(record.as_spectrum())
            if include_len_0:
                spectrum_inchikey_dict[inchikey] = spectrum_list
                # spectrum_smiles_dict[smiles].append(spectrum_list)
            elif include_len_0 is False:
                if len(spectrum_list):
                    spectrum_inchikey_dict[inchikey] = spectrum_list
                    # spectrum_smiles_dict[smiles].append(spectrum_list)
        return spectrum_inchikey_dict

    def get(self, spec_id):
        query = self.session.query(SpecInfo).filter(SpecInfo.spectrum_id.like(spec_id))
        return query.first().as_spectrum()

    def get_stream(self):
        query = self.session.query(SpecInfo)
        for result in query.all():
            yield result.as_spectrum()

    def get_rand(self, num):
        query = self.session.query(SpecInfo).order_by(func.random()).limit(num)
        for result in query.all():
            yield result.as_spectrum()

    def get_spectrum_by_id(self, spectrum_id):
        query = self.session.query(SpecInfo).filter(SpecInfo.spectrum_id == spectrum_id)
        result = query.first()
        if result:
            return result.as_spectrum()
        else:
            return None

    def get_spectrum_list_by_smiles(
        self, tags, smiles, instrument_type="%", precursor_type="%"
    ):
        query = self.session.query(SpecInfo).filter(
            SpecInfo.smiles.like(smiles),
            SpecInfo.instrument_type.like(instrument_type),
            SpecInfo.precursor_type.like(precursor_type),
            SpecInfo.tags.like("".join(["%", tags, "%"])),
        )
        spectrum_list = []
        for result in query.all():
            result_spectrum = result.as_spectrum()
            if len(result_spectrum.peaklist) > 4:
                spectrum_list.append(result.as_spectrum())
        return spectrum_list

    def get_spectrum_dict_by_smiles(
        self,
        tags,
        smiles_list,
        instrument_type,
        precursor_type,
        cut_smiles_by_empty_peaklist=False,
    ):
        """
        {smiles_0:[spectrum_smiles_0], smiles_1:[spectrum_smiles_1], ... ,smiles_n:[spectrum_smiles_n]}
        :param smiles_list:
        :param instrument_type:
        :param precursor_type:
        :return:
        """
        spectrum_smiles_dict = dict()
        for smiles in smiles_list:
            # print("get public spectrum dict one smiles!")
            spectrum_list = []
            query = self.session.query(SpecInfo).filter(
                SpecInfo.smiles == smiles,
                SpecInfo.instrument_type == instrument_type,
                SpecInfo.precursor_type == precursor_type,
                SpecInfo.tags == tags,
            )
            for result in query.all():
                result_spectrum = result.as_spectrum()
                if len(result_spectrum.peaklist) > 4:
                    spectrum_list.append(result.as_spectrum())
            if cut_smiles_by_empty_peaklist:
                if len(spectrum_list):
                    spectrum_smiles_dict[smiles] = spectrum_list
            elif cut_smiles_by_empty_peaklist is False:
                spectrum_smiles_dict[smiles] = spectrum_list
        return spectrum_smiles_dict

    def get_spectrum_list_by(self, smiles, instrutype, precutype):
        query = (
            self.session.query(SpecInfo)
            .filter(SpecInfo.smiles == smiles)
            .filter(SpecInfo.instrument_type == instrutype)
            .filter(SpecInfo.precursor_type == precutype)
        )
        spectrum_list = []
        for result in query.all():
            spectrum_list.append(result.as_spectrum())
        return spectrum_list

    def get_grouped_spectrum_by_smiles_stream(
        self, smiles, instrument_type, precursor_type, tags
    ):
        query = (
            self.session.query(
                SpecInfo.smiles,
                SpecInfo.instrument_type,
                SpecInfo.precursor_type,
                SpecInfo.author,
                SpecInfo.num_decimal_places,
                func.group_concat(SpecInfo.spectrum_id).label("group_id"),
            )
            .filter(
                SpecInfo.smiles == smiles,
                SpecInfo.instrument_type == instrument_type,
                SpecInfo.precursor_type == precursor_type,
                SpecInfo.tags == tags,
            )
            .group_by(
                SpecInfo.smiles,
                SpecInfo.instrument_type,
                SpecInfo.precursor_type,
                SpecInfo.author,
                SpecInfo.num_decimal_places,
            )
        )

        for result in query.all():
            spectrum_list = []
            for spectrum_id in result.group_id.split(","):
                try:
                    spectrum_list.append(self.get_spectrum_by_id(spectrum_id))
                except SyntaxError as se:
                    print(
                        "SyntaxError occurred while processing {!r}".format(spectrum_id)
                    )
                except ValueError as ve:
                    print(
                        "ValueError occurred while processing {!r}".format(spectrum_id)
                    )
            if len(spectrum_list) > 0:
                yield spectrum_list

    def get_grouped_spectrum_by_inchikey_stream(
        self, inchikey, instrument_type, precursor_type, tags
    ):
        query = (
            self.session.query(
                SpecInfo.inchikey,
                SpecInfo.instrument_type,
                SpecInfo.precursor_type,
                SpecInfo.author,
                SpecInfo.num_decimal_places,
                func.group_concat(SpecInfo.spectrum_id).label("group_id"),
            )
            .filter(
                SpecInfo.inchikey == inchikey,
                SpecInfo.instrument_type == instrument_type,
                SpecInfo.precursor_type == precursor_type,
                SpecInfo.tags == tags,
            )
            .group_by(
                SpecInfo.inchikey,
                SpecInfo.instrument_type,
                SpecInfo.precursor_type,
                SpecInfo.author,
                SpecInfo.num_decimal_places,
            )
        )

        for result in query.all():
            spectrum_list = []
            for spectrum_id in result.group_id.split(","):
                try:
                    spectrum_list.append(self.get_spectrum_by_id(spectrum_id))
                except SyntaxError as se:
                    print(
                        "SyntaxError occurred while processing {!r}".format(spectrum_id)
                    )
                except ValueError as ve:
                    print(
                        "ValueError occurred while processing {!r}".format(spectrum_id)
                    )
            if len(spectrum_list) > 0:
                yield spectrum_list

    def filtered_spectrum_group_stream(self, spec_filter=None):
        if not spec_filter:

            def default_filter(arg):
                return True

            spec_filter = default_filter

        query = self.session.query(
            SpecInfo.smiles,
            SpecInfo.instrument_type,
            SpecInfo.precursor_type,
            SpecInfo.author,
            SpecInfo.num_decimal_places,
            func.group_concat(SpecInfo.spectrum_id).label("group_id"),
        ).group_by(
            SpecInfo.smiles,
            SpecInfo.instrument_type,
            SpecInfo.precursor_type,
            SpecInfo.author,
            SpecInfo.num_decimal_places,
        )
        for result in query.all():
            spectrum_list = []
            for spectrum_id in result.group_id.split(","):
                try:
                    spec = self.get_spectrum_by_id(spectrum_id)
                except SyntaxError:
                    print(
                        "SyntaxError occurred while processing {!r}".format(spectrum_id)
                    )
                else:
                    if spec_filter(spec):
                        spectrum_list.append(spec)
            yield spectrum_list

    def get_grouped_count_stream(self, minimum_count=0):
        from sqlalchemy import distinct, desc

        query = (
            self.session.query(
                SpecInfo.smiles,
                SpecInfo.instrument_type,
                SpecInfo.precursor_type,
                func.count(distinct(SpecInfo.author)).label("C"),
            )
            .group_by(
                SpecInfo.smiles, SpecInfo.instrument_type, SpecInfo.precursor_type
            )
            .having(func.count(distinct(SpecInfo.author)) >= minimum_count)
            .order_by(desc("C"))
        )

        for result in query.all():
            result_dict = {
                "smiles": result.smiles,
                "instrument_type": result.instrument_type,
                "precursor_type": result.precursor_type,
                "count": result.C,
            }
            yield result_dict

    def count(self):
        return self.session.query(SpecInfo.spectrum_id).count()

    def truncate(self):
        self.session.execute("DELETE FROM info")
        self.session.commit()

    def close(self):
        self.session.close()
