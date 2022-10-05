import numpy as np
from sklearn.linear_model import LinearRegression
from pyopenms import MSExperiment, MzMLFile, Precursor, MSSpectrum
from data.chemical import Chemical


class Spectrum:
    __slots__ = [
        "_spec_id",
        "_chem_id",
        "_snapeaks_id",
        "_year",
        "_batch",
        "_batch_group",
        "_instrument_type",
        "_prec_type",
        "_prec_mz",
        "_collision_energy",
        "_peaks",
        "_noise_mean",
        "_noise_median",
        "_noise_mad",
        "_noise_dynamic",
        "_smiles",
        "_inchikey",
        "_metfrag_valid",
    ]

    def __init__(self):
        self._spec_id = None
        self._chem_id = None
        self._snapeaks_id = None
        self._prec_type = None
        self._prec_mz = None
        self._instrument_type = None
        self._collision_energy = None
        self._peaks = None
        self._noise_mean = None
        self._noise_median = None
        self._noise_mad = None
        self._noise_dynamic = None
        self._smiles = None
        self._inchikey = None
        self._metfrag_valid = None

    def __repr__(self):
        return (
            "<spectrum(spec_id={spec_id}, "
            "chem_id={chem_id}, "
            "snapeaks_id={snapeaks_id}, "
            "instrument_type={instrument_type}, "
            "prec_type={prec_type}, "
            "prec_mz={prec_mz}, "
            "collision_energy={collision_energy}, "
            "estimate_noise:"
            "mean={estimate_noise_mean},"
            "median={estimate_noise_median},"
            "mad={estimate_noise_mad},"
            "dynamic={estimate_noise_dynamic}, "
            "smiles={smiles}, "
            "inchikey={inchikey},"
            "metfrag_valid={metfrag_valid},"
            "peaks={peaks}".format(
                spec_id=self._spec_id,
                chem_id=self._chem_id,
                snapeaks_id=self._snapeaks_id,
                instrument_type=self._instrument_type,
                prec_type=self._prec_type,
                prec_mz=self._prec_mz,
                collision_energy=self._collision_energy,
                estimate_noise_mean=self._noise_mean,
                estimate_noise_median=self._noise_median,
                estimate_noise_mad=self._noise_mad,
                estimate_noise_dynamic=self._noise_dynamic,
                smiles=self._smiles,
                inchikey=self._inchikey,
                metfrag_valid=self._metfrag_valid,
                peaks=self._peaks,
            )
        )

    def info(self):
        return (
            "<spectrum(spec_id={spec_id}, "
            "chem_id={chem_id}, "
            "snapeaks_id={snapeaks_id}, "
            "instrument_type={instrument_type}, "
            "prec_type={prec_type}, "
            "prec_mz={prec_mz}, "
            "collision_energy={collision_energy}, "
            "estimate_noise:mean={noise_mean},"
            "median={noise_median},"
            "mad={noise_mad}, "
            "dynamic={noise_dynamic},"
            "smiles={smiles},"
            "inchikey={inchikey}, metfrag_valid={metfrag_valid}".format(
                spec_id=self._spec_id,
                chem_id=self._chem_id,
                snapeaks_id=self._snapeaks_id,
                instrument_type=self._instrument_type,
                prec_type=self._prec_type,
                prec_mz=self._prec_mz,
                collision_energy=self._collision_energy,
                noise_mean=self._noise_mean,
                noise_median=self._noise_median,
                noise_mad=self._noise_mad,
                noise_dynamic=self._noise_dynamic,
                smiles=self._smiles,
                inchikey=self._inchikey,
                metfrag_valid=self._metfrag_valid,
            )
        )

    @property
    def spec_id(self):
        return self._spec_id

    @spec_id.setter
    def spec_id(self, spec_id):
        self._spec_id = spec_id

    @property
    def chem_id(self):
        return self._chem_id

    @chem_id.setter
    def chem_id(self, chem_id):
        self._chem_id = chem_id

    @property
    def snapeaks_id(self):
        return self._snapeaks_id

    @snapeaks_id.setter
    def snapeaks_id(self, snapeaks_id):
        self._snapeaks_id = snapeaks_id

    @property
    def instrument_type(self):
        return self._instrument_type

    @instrument_type.setter
    def instrument_type(self, instrument_type):
        self._instrument_type = instrument_type

    @property
    def prec_type(self):
        return self._prec_type

    @prec_type.setter
    def prec_type(self, prec_type):
        self._prec_type = prec_type

    @property
    def prec_mz(self):
        return self._prec_mz

    @prec_mz.setter
    def prec_mz(self, prec_mz):
        self._prec_mz = prec_mz

    @property
    def collision_energy(self):
        return self._collision_energy

    @collision_energy.setter
    def collision_energy(self, collision_energy):
        self._collision_energy = collision_energy

    @property
    def peaks(self):
        return self._peaks

    @peaks.setter
    def peaks(self, peaks):
        self._peaks = peaks

    @property
    def noise_mean(self):
        return self._noise_mean

    @noise_mean.setter
    def noise_mean(self, noise_mean):
        self._noise_mean = noise_mean

    @property
    def noise_median(self):
        return self._noise_median

    @noise_median.setter
    def noise_median(self, noise_median):
        self._noise_median = noise_median

    @property
    def noise_mad(self):
        return self._noise_mad

    @noise_mad.setter
    def noise_mad(self, noise_mad):
        self._noise_mad = noise_mad

    @property
    def noise_dynamic(self):
        return self._noise_dynamic

    @noise_dynamic.setter
    def noise_dynamic(self, noise_dynamic):
        self._noise_dynamic = noise_dynamic

    @property
    def smiles(self):
        return self._smiles

    @smiles.setter
    def smiles(self, smiles):
        self._smiles = smiles

    @property
    def inchikey(self):
        return self._inchikey

    @inchikey.setter
    def inchikey(self, inchikey):
        self._inchikey = inchikey

    @property
    def metfrag_valid(self):
        return self._metfrag_valid

    @metfrag_valid.setter
    def metfrag_valid(self, metfrag_valid):
        self._metfrag_valid = metfrag_valid

    def estimate_noise_dynamic(self):
        intensity_np_list = self._peaks[:, 1]
        np_id = np.array([[i] for i in range(1, len(intensity_np_list) + 1)])
        intensity_np_list = np.sort(intensity_np_list)
        np_intensity_reshape = intensity_np_list.reshape(-1, 1)
        estimate_dynamic_noise = np.max(intensity_np_list)  # init value
        snr = 2.0
        estimate_dynamic_noise = np.max(intensity_np_list)
        for n in range(3, len(np_intensity_reshape)):
            np_id_set = np_id[:n]
            np_intensity_reshape_set = np_intensity_reshape[:n]  # extract array 0 ~ n-1
            model = LinearRegression()
            model.fit(X=np_id_set, y=np_intensity_reshape_set)
            # print(np_intensity_reshape[n][0], model.predict([[n]])[0][0],
            #       np_intensity_reshape[n][0] / model.predict([[n]])[0][0],
            #       (np_intensity_reshape[n][0] / model.predict([[n]])[0][0]) > snr)
            if np_intensity_reshape[n][0] / model.predict([[n]])[0][0] > snr:
                estimate_dynamic_noise = np_intensity_reshape[n - 1][0]
                break
        return estimate_dynamic_noise

    def peaks_cut_noise(self, mod):
        temp_spec = self.peaks.copy()
        if mod == "original":
            cut_intensity = np.max(temp_spec[:, 1]) * 0.05
            # temp_spec = temp_spec[temp_spec[:, 1] > cut_intensity]  # cut off 5%
            temp_spec[:, 1] = temp_spec[:, 1] - cut_intensity
        elif mod == "mean":
            # temp_spec = temp_spec[temp_spec[:, 1] > self.noise_mean, :]
            temp_spec[:, 1] = temp_spec[:, 1] - self.noise_mean
        elif mod == "mad":
            # temp_spec = temp_spec[temp_spec[:, 1] > self.noise_mad, :]
            temp_spec[:, 1] = temp_spec[:, 1] - self.noise_mad
        elif mod == "median":
            # temp_spec = temp_spec[temp_spec[:, 1] > self.noise_median, :]
            temp_spec[:, 1] = temp_spec[:, 1] - self.noise_median
        elif mod == "dynamic":
            # temp_spec = temp_spec[temp_spec[:, 1] > self.noise_dynamic, :]
            temp_spec[:, 1] = temp_spec[:, 1] - self.noise_dynamic

        temp_spec = temp_spec[temp_spec[:, 1] > 0.0, :]
        if len(temp_spec) > 0:
            max_intensity = np.max(temp_spec[:, 1])
            temp_spec[:, 1] = temp_spec[:, 1] / max_intensity * 1000.0
            return temp_spec
        else:
            return np.array([[]])

    def dynamic_peaks(self):
        temp_modify_spec = self.peaks
        temp_spec = self.peaks
        max_intensity = np.max(temp_spec[:, 1])
        if len(temp_spec) > 2:
            smallest_intensity = np.partition(temp_spec[:, 1], 2)[0]
            second_smallest_intensity = np.partition(temp_spec[:, 1], 2)[1]
            if second_smallest_intensity - smallest_intensity < 0.05 * max_intensity:
                temp_spec[:, 1] = temp_spec[:, 1] - second_smallest_intensity
            else:
                temp_spec[:, 1] = temp_spec[:, 1] - smallest_intensity
            temp_spec = temp_spec[temp_spec[:, 1] > 0.0, :]
            if len(temp_spec):
                dynamic_noise_level = self.calculate_noise_dynamic(temp_spec)
                temp_modify_spec = temp_modify_spec[
                    temp_modify_spec[:, 1] > dynamic_noise_level, :
                ]
        else:
            temp_modify_spec = temp_modify_spec[
                temp_modify_spec[:, 1] > self.noise_dynamic, :
            ]
        temp_modify_spec = temp_modify_spec[temp_modify_spec[:, 1] > 0.0, :]
        if len(temp_modify_spec):
            max_intensity = np.max(temp_modify_spec[:, 1])
            temp_modify_spec[:, 1] = temp_modify_spec[:, 1] / max_intensity * 1000.0
            return temp_modify_spec
        else:
            return np.array([[]])

    @staticmethod
    def calculate_noise_dynamic(peaks):
        intensity_np_list = peaks[:, 1]
        np_id = np.array([[i] for i in range(1, len(intensity_np_list) + 1)])
        intensity_np_list = np.sort(intensity_np_list)
        np_intensity_reshape = intensity_np_list.reshape(-1, 1)
        estimate_dynamic_noise = np.max(intensity_np_list)  # init value
        snr = 1.5
        estimate_dynamic_noise = np.max(intensity_np_list)
        for n in range(3, len(np_intensity_reshape)):
            np_id_set = np_id[:n]
            np_intensity_reshape_set = np_intensity_reshape[:n]  # extract array 0 ~ n-1
            model = LinearRegression()
            model.fit(X=np_id_set, y=np_intensity_reshape_set)
            # print(np_intensity_reshape[n][0], model.predict([[n]])[0][0],
            #       np_intensity_reshape[n][0] / model.predict([[n]])[0][0],
            #       (np_intensity_reshape[n][0] / model.predict([[n]])[0][0]) > snr)
            if np_intensity_reshape[n][0] / model.predict([[n]])[0][0] > snr:
                estimate_dynamic_noise = np_intensity_reshape[n - 1][0]
                break
        return estimate_dynamic_noise

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__slots__ == other.__slots__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    @staticmethod
    def build_spec(
        spec_id: str,
        precursor: Precursor,
        precursor_type,
        spectrum: MSSpectrum,
        compound: Chemical,
        instrument_type: str
    ):
        spec = Spectrum()
        spec.spec_id = spec_id
        spec.c_id = compound.c_id
        spec.bmdms_id = compound.bmdms_id
        spec.native_id_in_mzml = spectrum.getNativeID()
        spec.instrument_type = instrument_type
        spec.retention_time = spectrum.getRT()
        spec.prec_type = precursor_type
        spec.prec_mz = precursor.getMZ()
        spec.collision_energy = precursor.getMetaValue("collision energy")
        spec.ms_level = spectrum.getMSLevel()
        spec.tic = spectrum.getTIC()
        peak_list = []
        for p in spectrum:
            peak_list.append([p.getMZ(), p.getIntensity()])
        spec.peaks = np.array(peak_list)
        spec.peaks = spec.peaks[spec.peaks[:, 1] != 0]
        spec.noise_dynamic = (
            spec.estimate_noise_dynamic()
        )  # 노이즈 제거가 여기에 들어가는 게 맞을지. MS3 peak 추출 후에 noise 제거가 되어야할 수 있음.
        return spec