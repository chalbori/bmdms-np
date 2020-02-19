import os
import sys, inspect

CURRENT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, PARENT_DIR)
print("sys.path: {}".format(sys.path))
import tools.tool_file as tool_file
from database import info_chem_mysql, info_spectrum_mysql
import data.spectrum as spectrum
from database.accounts import Database
import numpy as np
import time

smi_list = []
db_uri = 'mysql+pymysql://{}:{}@{}/{}'.format(Database.user_id, Database.user_password, Database.host,
                                              Database.database)
repo_spectrum = info_spectrum_mysql.SpectrumRepository(db_uri=db_uri)
repo_chemical = info_chem_mysql.CompoundRepository(db_uri=db_uri)

if __name__ == '__main__':
    from spectrum_utils.spectrum import MsmsSpectrum
    import spectrum_utils.plot as plot
    import matplotlib.pyplot as plt

    spectrum_stream = repo_spectrum.get_stream(include_null=False, limit=100)
    count = 0
    for spectrum in spectrum_stream:
        # print(spectrum.spec_id)
        peaks_nr = spectrum.peaks_cut_noise(mod='dynamic')
        if peaks_nr.shape[1]:
            temp = MsmsSpectrum(identifier=spectrum.spec_id, precursor_mz=spectrum.prec_mz, precursor_charge=1,
                                mz=peaks_nr[:, 0],
                                intensity=peaks_nr[:, 1])

            plt.figure(figsize=(6, 4), dpi=300)
            ax = plot.spectrum(temp)

            output_folder = 'out/spectra/{}_{}/{}/{}/'.format(spectrum.snapeaks_id, spectrum.inchikey,
                                                              spectrum.instrument_type, spectrum.prec_type)
            tool_file.mkdir_folder(output_folder)

            # i = 0
            output_file = output_folder + '{}eV_specid_{}.png'.format(spectrum.collision_energy, spectrum.spec_id)
            # while os.path.exists(output_file):
            #     i += 1
            #     output_file = output_folder + '{}eV_specid_{}.png'.format(spectrum.collision_energy, spectrum.spec_id)

            plt.savefig(output_file)
            plt.close()
        count += 1
        if count % 1000 == 0:
            print('count : {}'.format(count), time.strftime('%c', time.localtime(time.time())))


    