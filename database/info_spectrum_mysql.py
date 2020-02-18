import io
import sqlite3
import sqlalchemy
from sqlalchemy import Column, String, Integer, Float, LargeBinary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import func
from sqlalchemy.dialects import mysql
import numpy as np
from data.spectrum import Spectrum

Base = declarative_base()


class SpectrumInfo(Base):
    __tablename__ = "info_spectrum"

    spec_id = Column(Integer, autoincrement=True, primary_key=True)
    chem_id = Column(Integer)
    snapeaks_id = Column(String(20))
    prec_type = Column(String(15))
    prec_mz = Column(Float)
    collision_energy = Column(Integer)
    peaks = Column(LargeBinary(length=(2 ** 32) - 1))
    instrument_type = Column(String(15))
    estimate_noise_mean = Column(Float)
    estimate_noise_median = Column(Float)
    estimate_noise_mad = Column(Float)
    estimate_noise_dynamic = Column(Float)
    smiles = Column(String(5000))
    inchikey = Column(mysql.VARCHAR(30))
    metfrag_valid = Column(mysql.BOOLEAN)

    @staticmethod
    def from_spectrum(spectrum):
        spec_info = SpectrumInfo()
        spec_info.spec_id = spectrum.spec_id
        spec_info.chem_id = spectrum.chem_id
        spec_info.snapeaks_id = spectrum.snapeaks_id
        spec_info.prec_type = spectrum.prec_type
        spec_info.prec_mz = spectrum.prec_mz
        spec_info.collision_energy = spectrum.collision_energy
        spec_info.peaks = adapt_array(spectrum.peaks)
        spec_info.instrument_type = spectrum.instrument_type
        spec_info.estimate_noise_mean = spectrum.noise_mean
        spec_info.estimate_noise_median = spectrum.noise_median
        spec_info.estimate_noise_mad = spectrum.noise_mad
        spec_info.estimate_noise_dynamic = spectrum.noise_dynamic
        spec_info.smiles = spectrum.smiles
        spec_info.inchikey = spectrum.inchikey
        spec_info.metfrag_valid = spectrum.metfrag_valid
        return spec_info

    def as_spectrum(self):
        _spectrum = Spectrum()
        _spectrum.spec_id = self.spec_id
        _spectrum.chem_id = self.chem_id
        _spectrum.snapeaks_id = self.snapeaks_id
        _spectrum.prec_type = self.prec_type
        _spectrum.prec_mz = self.prec_mz
        _spectrum.instrument_type = self.instrument_type
        _spectrum.collision_energy = self.collision_energy
        _spectrum.peaks = convert_array(self.peaks)
        _spectrum.noise_mean = self.estimate_noise_mean
        _spectrum.noise_median = self.estimate_noise_median
        _spectrum.noise_mad = self.estimate_noise_mad
        _spectrum.noise_dynamic = self.estimate_noise_dynamic
        _spectrum.smiles = self.smiles
        _spectrum.inchikey = self.inchikey
        _spectrum.metfrag_valid = self.metfrag_valid
        return _spectrum

    def __repr__(self):
        return (
            "<spectrum_info(spec_id={spec_id}, "
            "chem_id={chem_id}, "
            "snapeaks_id={snapeaks_id}, "
            "prec_type={prec_type}, "
            "prec_mz={prec_mz}, "
            "collision_energy={collision_energy}, "
            "instrument_type={instrument_type}, "
            "estimate_noise_mean={estimate_noise_mean}, "
            "estimate_noise_median={estimate_noise_median}, "
            "estimate_noise_mad={estimate_noise_mad}, "
            "estimate_noise_dynamic={estimate_noise_dynamic}, "
            "smiles={smiles},"
            "inchikey={inchikey}, metfrag_valid={metfrag_valid}".format(
                spec_id=self.spec_id,
                chem_id=self.chem_id,
                snapeaks_id=self.snapeaks_id,
                prec_type=self.prec_type,
                prec_mz=self.prec_mz,
                collision_energy=self.collision_energy,
                instrument_type=self.instrument_type,
                estimate_noise_mean=self.estimate_noise_mean,
                estimate_noise_median=self.estimate_noise_median,
                estimate_noise_mad=self.estimate_noise_mad,
                estimate_noise_dynamic=self.estimate_noise_dynamic,
                smiles=self.smiles,
                inchikey=self.inchikey,
                metfrag_valid=self.metfrag_valid,
            )
        )


class SpectrumRepository:
    def __init__(self, db_uri):
        engine = sqlalchemy.create_engine(db_uri, echo=False)
        Base.metadata.create_all(engine)
        Session = sessionmaker()
        Session.configure(bind=engine)
        self.session = Session()

    def get(self, spec_id):
        query = self.session.query(SpectrumInfo).filter(SpectrumInfo.spec_id == spec_id)
        return query.first().as_spectrum()

    def getrand(self, spec_id="%", sample_num=100):
        query = (
            self.session.query(SpectrumInfo)
            .filter(SpectrumInfo.spec_id.like(spec_id))
            .order_by(func.random())
            .limit(sample_num)
        )
        for record in query.all():
            yield record.as_spectrum()

    def get_rand(self, spec_id="%", sample_num=100, metfrag_valid="%"):
        query = (
            self.session.query(SpectrumInfo)
            .filter(SpectrumInfo.spec_id.like(spec_id), SpectrumInfo.metfrag_valid.like(metfrag_valid))
            .order_by(func.random())
            .limit(sample_num)
        )
        for record in query.all():
            yield record.as_spectrum()

    def gettest(
        self,
        inchikeylist,
        instrument_type,
        precursor_type,
        collision_energy,
        sample_num,
    ):

        spectrum_inchikey_dict = dict()
        for inchikey in inchikeylist:
            temp = []
            query = (
                self.session.query(SpectrumInfo)
                .filter(
                    SpectrumInfo.inchikey == inchikey,
                    SpectrumInfo.instrument_type == instrument_type,
                    SpectrumInfo.prec_type == precursor_type,
                    SpectrumInfo.collision_energy == collision_energy,
                )
                .order_by(func.rand())
                .limit(int(sample_num / 5))
            )
            for record in query.all():
                temp.append(record.as_spectrum())
            spectrum_inchikey_dict[inchikey] = temp
        return spectrum_inchikey_dict

    def get_id_peaks(self, spec_id="%"):
        query = self.session.query(SpectrumInfo.spec_id, SpectrumInfo.peaks).filter(
            SpectrumInfo.spec_id == spec_id
        )
        return query.first().spec_id, query.first().peaks

    def get_stream(self, include_null=False, limit=10000000):
        if include_null is True:
            query = (
                self.session.query(SpectrumInfo)
                .order_by(SpectrumInfo.spec_id)
                .limit(limit=limit)
            )
        else:
            query = (
                self.session.query(SpectrumInfo)
                .filter(SpectrumInfo.snapeaks_id.isnot(None))
                .order_by(SpectrumInfo.spec_id)
                .limit(limit=limit)
            )
        for result in query.all():
            yield result.as_spectrum()

    def get_stream_range(self, include_null=False, range=[0, 1000000]):
        if include_null is True:
            query = (
                self.session.query(SpectrumInfo)
                .filter(
                    SpectrumInfo.spec_id >= range[0], SpectrumInfo.spec_id <= range[1]
                )
                .order_by(SpectrumInfo.spec_id)
            )
        else:
            query = (
                self.session.query(SpectrumInfo)
                .filter(
                    SpectrumInfo.snapeaks_id.isnot(None),
                    SpectrumInfo.spec_id >= range[0],
                    SpectrumInfo.spec_id <= range[1],
                )
                .order_by(SpectrumInfo.spec_id)
            )
        for result in query.all():
            yield result.as_spectrum()

    def get_smiles_list(self, instrument_type="%", precursor_type="%"):
        smiles_list = []
        query = (
            self.session.query(SpectrumInfo.smiles)
            .distinct(SpectrumInfo.smiles)
            .filter(
                SpectrumInfo.instrument_type.like(instrument_type),
                SpectrumInfo.prec_type.like(precursor_type),
            )
        )
        for record in query.all():
            smiles_list.append(record.smiles)
        return smiles_list

    def get_smiles_list_metfrag_valid(self, instrument_type="%", precursor_type="%"):
        smiles_list = []
        query = (
            self.session.query(SpectrumInfo.smiles)
            .distinct(SpectrumInfo.smiles)
            .filter(
                SpectrumInfo.instrument_type.like(instrument_type),
                SpectrumInfo.prec_type.like(precursor_type),
                SpectrumInfo.metfrag_valid == 1,
            )
        )
        for record in query.all():
            smiles_list.append(record.smiles)
        return smiles_list

    def get_inchikey_list_metfrag_valid(self, instrument_type="%", precursor_type="%"):
        inchikey_list = []
        query = (
            self.session.query(SpectrumInfo.inchikey)
            .distinct(SpectrumInfo.inchikey)
            .filter(
                SpectrumInfo.instrument_type.like(instrument_type),
                SpectrumInfo.prec_type.like(precursor_type),
                SpectrumInfo.metfrag_valid == 1,
            )
        )
        for record in query.all():
            inchikey_list.append(record.inchikey)
        return inchikey_list

    def get_cid_list(self, instrument_type="%", precursor_type="%"):
        smiles_list = []
        query = (
            self.session.query(SpectrumInfo.snapeaks_id)
            .distinct(SpectrumInfo.snapeaks_id)
            .filter(
                SpectrumInfo.instrument_type.like(instrument_type),
                SpectrumInfo.prec_type.like(precursor_type),
                SpectrumInfo.snapeaks_id.isnot(None),
            )
        )
        for record in query.all():
            smiles_list.append(record.snapeaks_id)
        return smiles_list

    def get_spectrum_dict_by_smiles(
        self, smiles_list, instrument_type="%", precursor_type="%", include_len_0=True
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
            # print("get spectrum dict one smiles!")
            spectrum_list = []
            query = self.session.query(SpectrumInfo).filter(
                SpectrumInfo.smiles == smiles,
                SpectrumInfo.instrument_type == instrument_type,
                SpectrumInfo.prec_type == precursor_type,
            )
            for result in query.all():
                result_spectrum = result.as_spectrum()
                if len(result_spectrum.peaks) > 4:
                    spectrum_list.append(result.as_spectrum())
            if include_len_0:
                spectrum_smiles_dict[smiles] = spectrum_list
            elif include_len_0 is False:
                if len(spectrum_list):
                    spectrum_smiles_dict[smiles] = spectrum_list
        return spectrum_smiles_dict

    def get_spectrum_dict_by_smiles_choose(
        self,
        smiles_list,
        instrument_type="%",
        precursor_type="%",
        include_len_0=True,
        sample_num=20,
    ):
        """
        {smiles_0:[spectrum_smiles_0], smiles_1:[spectrum_smiles_1], ... ,smiles_n:[spectrum_smiles_n]}
        :param smiles_list:
        :param instrument_type:
        :param precursor_type:
        :return:
        """
        collision_energy_list = [10, 20, 40, 60, 80]
        spectrum_smiles_dict = dict()
        for smiles in smiles_list:
            # print("get spectrum dict one smiles!")
            spectrum_list = []
            for collision_energy in collision_energy_list:
                query = (
                    self.session.query(SpectrumInfo)
                    .filter(
                        SpectrumInfo.smiles == smiles,
                        SpectrumInfo.instrument_type == instrument_type,
                        SpectrumInfo.prec_type == precursor_type,
                        SpectrumInfo.collision_energy == collision_energy,
                    )
                    .order_by(func.random())
                    .limit(int(sample_num / len(collision_energy_list)))
                )
                for result in query.all():
                    result_spectrum = result.as_spectrum()
                    if len(result_spectrum.peaks) > 4:
                        spectrum_list.append(result.as_spectrum())
            if include_len_0:
                spectrum_smiles_dict[smiles] = spectrum_list
                # spectrum_smiles_dict[smiles].append(spectrum_list)
            elif include_len_0 is False:
                if len(spectrum_list):
                    spectrum_smiles_dict[smiles] = spectrum_list
                    # spectrum_smiles_dict[smiles].append(spectrum_list)
        return spectrum_smiles_dict

    def get_spectrum_dict_by_smiles_choose_metfrag_valid(
        self,
        smiles_list,
        instrument_type="%",
        precursor_type="%",
        include_len_0=True,
        sample_num=20,
    ):
        """
        {smiles_0:[spectrum_smiles_0], smiles_1:[spectrum_smiles_1], ... ,smiles_n:[spectrum_smiles_n]}
        :param smiles_list:
        :param instrument_type:
        :param precursor_type:
        :return:
        """
        collision_energy_list = [10, 20, 40, 60, 80]
        spectrum_smiles_dict = dict()
        for smiles in smiles_list:
            # print("get spectrum dict one smiles!")
            spectrum_list = []
            for collision_energy in collision_energy_list:
                query = (
                    self.session.query(SpectrumInfo)
                    .filter(
                        SpectrumInfo.smiles == smiles,
                        SpectrumInfo.instrument_type == instrument_type,
                        SpectrumInfo.prec_type == precursor_type,
                        SpectrumInfo.collision_energy == collision_energy,
                        SpectrumInfo.metfrag_valid == 1,
                    )
                    .order_by(func.random())
                    .limit(int(sample_num / len(collision_energy_list)))
                )
                for result in query.all():
                    result_spectrum = result.as_spectrum()
                    if len(result_spectrum.peaks) > 4:
                        spectrum_list.append(result.as_spectrum())
            if include_len_0:
                spectrum_smiles_dict[smiles] = spectrum_list
                # spectrum_smiles_dict[smiles].append(spectrum_list)
            elif include_len_0 is False:
                if len(spectrum_list):
                    spectrum_smiles_dict[smiles] = spectrum_list
                    # spectrum_smiles_dict[smiles].append(spectrum_list)
        return spectrum_smiles_dict

    def get_spectrum_dict_by_inchikey_choose(
        self,
        inchikey_list,
        instrument_type="%",
        precursor_type="%",
        include_len_0=True,
        sample_num=20,
    ):
        collision_energy_list = [10, 20, 40, 60, 80]
        sampling_num_each = int(sample_num / len(collision_energy_list))
        spectrum_inchikey_dict = dict()
        for inchikey in inchikey_list:
            # print("get spectrum dict one smiles!")
            spectrum_list = []
            for collision_energy in collision_energy_list:
                # func.rand() : mysql or mariadb
                query = self.session.query(SpectrumInfo).filter(
                    SpectrumInfo.inchikey == inchikey,
                    SpectrumInfo.instrument_type == instrument_type,
                    SpectrumInfo.prec_type == precursor_type,
                    SpectrumInfo.collision_energy == collision_energy,
                )
                query = query.order_by(func.rand()).limit(sampling_num_each)
                for record in query.all():
                    result_spectrum = record.as_spectrum()
                    if len(result_spectrum.peaks) > 4:
                        spectrum_list.append(record.as_spectrum())
            if include_len_0:
                spectrum_inchikey_dict[inchikey] = spectrum_list
                # spectrum_smiles_dict[smiles].append(spectrum_list)
            elif include_len_0 is False:
                if len(spectrum_list) > 0:
                    spectrum_inchikey_dict[inchikey] = spectrum_list
                    # spectrum_smiles_dict[smiles].append(spectrum_list)
        return spectrum_inchikey_dict

    def get_spectrum_dict_by_inchikey_choose_metfrag_valid(
        self,
        inchikey_list,
        instrument_type="%",
        precursor_type="%",
        include_len_0=True,
        sample_num=20,
    ):
        """
        {inchikey_0:[spectrum_inchikey_0], inchikey_1:[spectrum_inchikey_1], ... ,inchikey_n:[spectrum_inchikey_n]}
        :param inchikey_list:
        :param instrument_type:
        :param precursor_type:
        :return:
        """
        collision_energy_list = [10, 20, 40, 60, 80]
        spectrum_inchikey_dict = dict()
        for inchikey in inchikey_list:
            # print("get spectrum dict one inchikey!")
            spectrum_list = []
            for collision_energy in collision_energy_list:
                query = (
                    self.session.query(SpectrumInfo)
                    .filter(
                        SpectrumInfo.inchikey == inchikey,
                        SpectrumInfo.instrument_type == instrument_type,
                        SpectrumInfo.prec_type == precursor_type,
                        SpectrumInfo.collision_energy == collision_energy,
                        SpectrumInfo.metfrag_valid == 1,
                    )
                    .order_by(func.random())
                    .limit(int(sample_num / len(collision_energy_list)))
                )
                for result in query.all():
                    result_spectrum = result.as_spectrum()
                    if len(result_spectrum.peaks) > 4:
                        spectrum_list.append(result.as_spectrum())
            if include_len_0:
                spectrum_inchikey_dict[inchikey] = spectrum_list
                # spectrum_inchikey_dict[inchikey].append(spectrum_list)
            elif include_len_0 is False:
                if len(spectrum_list):
                    spectrum_inchikey_dict[inchikey] = spectrum_list
                    # spectrum_inchikey_dict[inchikey].append(spectrum_list)
        return spectrum_inchikey_dict

    def get_spectrum_dict_by_cid_choose(
        self,
        cid_list,
        instrument_type="%",
        precursor_type="%",
        include_len_0=True,
        sample_num=20,
    ):
        """

        :param smiles_list:
        :param instrument_type:
        :param precursor_type:
        :return:
        """
        collision_energy_list = [10, 20, 40, 60, 80]
        spectrum_cid_dict = dict()
        for cid in cid_list:
            print("get spectrum id:{}!".format(cid))
            spectrum_list = []
            for collision_energy in collision_energy_list:
                query = (
                    self.session.query(SpectrumInfo)
                    .filter(
                        SpectrumInfo.snapeaks_id == cid,
                        SpectrumInfo.instrument_type == instrument_type,
                        SpectrumInfo.prec_type == precursor_type,
                        SpectrumInfo.collision_energy == collision_energy,
                    )
                    .order_by(func.random())
                    .limit(int(sample_num / len(collision_energy_list)))
                )
                for result in query.all():
                    result_spectrum = result.as_spectrum()
                    if len(result_spectrum.peaks) > 4:
                        spectrum_list.append(result.as_spectrum())
            if include_len_0:
                spectrum_cid_dict[cid] = spectrum_list
            elif include_len_0 is False:
                if len(spectrum_list):
                    spectrum_cid_dict[cid] = spectrum_list
        return spectrum_cid_dict

    def get_spectrum_dict_by_cid_choose_draw(
        self,
        cid_list,
        instrument_type="%",
        precursor_type="%",
        include_len_0=True,
        sample_num=20,
    ):
        """
        :param smiles_list:
        :param instrument_type:
        :param precursor_type:
        :return:
        """
        collision_energy_list = [10, 20, 40, 60, 80]
        spectrum_cid_dict = dict()
        for cid in cid_list:
            print("get spectrum id:{}!".format(cid))
            spectrum_list = []
            for collision_energy in collision_energy_list:
                query = (
                    self.session.query(SpectrumInfo)
                    .filter(
                        SpectrumInfo.snapeaks_id == cid,
                        SpectrumInfo.instrument_type == instrument_type,
                        SpectrumInfo.prec_type == precursor_type,
                        SpectrumInfo.collision_energy == collision_energy,
                    )
                    .order_by(func.random())
                    .limit(int(sample_num / len(collision_energy_list)))
                )
                for result in query.all():
                    spectrum_list.append(result.as_spectrum())
            if include_len_0:
                spectrum_cid_dict[cid] = spectrum_list
            elif include_len_0 is False:
                if len(spectrum_list):
                    spectrum_cid_dict[cid] = spectrum_list
        return spectrum_cid_dict

    def get_grouped_spectrum_by_inchikey_stream(self, noise_mod='dynamic', list_length_min=3):
        query = (
            self.session.query(
                SpectrumInfo.inchikey,
                SpectrumInfo.instrument_type,
                SpectrumInfo.prec_type,
                SpectrumInfo.collision_energy,
                func.group_concat(SpectrumInfo.spec_id).label("group_id"),
            )
            .filter(
                SpectrumInfo.collision_energy.in_([10, 20, 40, 60, 80])
            )
            .group_by(
                SpectrumInfo.inchikey,
                SpectrumInfo.instrument_type,
                SpectrumInfo.prec_type,
                SpectrumInfo.collision_energy
            )
        )
        for result in query.all():
            spectrum_list = []
            for spectrum_id in result.group_id.split(","):
                temp_spectrum = self.get(spec_id = spectrum_id)
                if len(temp_spectrum.peaks_cut_noise(noise_mod)) > 4:
                    spectrum_list.append(self.get(spec_id=spectrum_id))
            if len(spectrum_list) >= list_length_min:
                yield result.inchikey, result.instrument_type, result.prec_type, result.collision_energy, spectrum_list

    def add(self, spectrum_list):
        for spectrum in spectrum_list:
            if not issubclass(type(spectrum), Spectrum):
                raise TypeError(
                    "{!r} is an unsupported type\n"
                    "{!r} is expected.".format(type(spectrum), Spectrum)
                )
            else:
                pass

        wrapped_spectrum_list = []
        for spectrum in spectrum_list:
            wrapped_spectrum = SpectrumInfo.from_spectrum(spectrum=spectrum)
            wrapped_spectrum_list.append(wrapped_spectrum)
        self.session.bulk_save_objects(wrapped_spectrum_list)
        self.session.commit()
        # if not issubclass(type(spectrum), Spectrum):
        #     raise TypeError('{!r} is an unsupported type\n'
        #                     '{!r} is expected.'.format(type(spectrum), Spectrum))
        # else:
        #     pass
        #
        # wrapped_spectrum = SpectrumInfo.from_spectrum(spectrum=spectrum)
        #
        # if not self.get(spec_id=spectrum.spec_id):
        #     self.session.add(wrapped_spectrum)
        #     self.session.commit()
        # else:
        #     pass

    def update_dynamic_noise_new(self, spec_id, value):
        # for spec in self.get('1'):
        # for spec in self.get_stream(include_null=True):
        self.session.query(SpectrumInfo).filter(SpectrumInfo.spec_id == spec_id).update(
            {SpectrumInfo.estimate_noise_dynamic: value}
        )
        self.session.commit()

    def update_dynamic_noise(self, mappings):
        # for spec in self.get('1'):
        # for spec in self.get_stream(include_null=True):
        self.session.bulk_update_mappings(SpectrumInfo, mappings=mappings)
        self.session.flush()
        self.session.commit()

    def update(self, mappings):
        self.session.bulk_update_mappings(SpectrumInfo, mappings=mappings)
        self.session.flush()
        self.session.commit()

    def close(self):
        self.session.close()


def adapt_array(arr):
    """
    http://stackoverflow.com/a/31312102/190597 (SoulNibbler)
    """
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return sqlite3.Binary(out.read())


def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)
