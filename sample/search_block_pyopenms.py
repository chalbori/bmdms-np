import os
import sys, inspect

CURRENT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, PARENT_DIR)
# print("sys.path: {}".format(sys.path))
import datetime
import multiprocessing
from pyopenms import MSExperiment, MzMLFile
import numpy as np
from data import spectrum
from database import info_spectrum_mysql
from database.accounts import Database
from tools import tool_file, tool_calculate, tool_list, tool_mzml
from tqdm import tqdm


NUM_PROCESSOR = multiprocessing.cpu_count() - 1
MULTI_PROC_UNIT = NUM_PROCESSOR * 5000
TOP_K_NUM = 20
NOISE_MOD = "dynamic"
URI_BMDMS = "mysql+pymysql://{}:{}@{}/{}".format(
        Database.user_id, Database.user_password, Database.host, Database.database
    )

def main(query_path_dir, instrument_type, ion_mode, precursor_mz):
    """
    Basic example script to search query block against BMDMS-NP

    usage:

        (source code folder/sample):~$ python search_block.py <query_path_dir> <instrument_type> <ion_mode> <precursor_mz> 

    example:
        (source code folder/sample):~$ python search_block.py ~/dev_vsc/bmdms-np/test/search_query_block/No6_RNXZPKOEJUFJON-UHFFFAOYSA-N/Orbitrap "Orbitrap" "[M+H]+" 331.082

    parameters:
        <instrument_type> : "Orbitrap", "O", "QTOF", "Q-TOF", "Q"
        <ion_mode> : "[M+H]+", "H", "[M+Na]+", "Na"
        <precursor_mz> : float (see compounds.csv)
    """
    
    repo_spectrum = info_spectrum_mysql.SpectrumRepository(db_uri=URI_BMDMS)
    file_list = tool_file.list_dir_except_hidden(path=query_path_dir)
    print(file_list)
    file_list.sort()

    print("INSTRUMENT TYPE: {}".format(instrument_type))
    print("ION MODE: {}".format(ion_mode))
    print("PRECURSOR MZ: {}".format(precursor_mz))

    query_spec_list = []
    ms_exp_list = tool_mzml.read_mzml(input_folder=query_path_dir)
    
    print("Extract mzML in query file(s)")
    for exp in tqdm(ms_exp_list):
        for s in exp.getSpectra():
            if s.getMSLevel() == 2:
                precursor = s.getPrecursors()[0]
                if abs(precursor.getMZ() - precursor_mz) < 0.05:
                
                    spec = spectrum.Spectrum()
                    spec.prec_mz = precursor.getMZ()
                    # prec_list.append(spec.prec_mz)
                    spec.instrument_type = 'Orbitrap'
                    # spec.collision_energy = int(collision_energy)
                    peak_list = []
                    for p in s:
                        peak_list.append([p.getMZ(), p.getIntensity()])
                    spec.peaks = np.array(peak_list)
                    spec.peaks = spec.peaks[spec.peaks[:, 1] != 0]
                    spec.peaks = spec.dynamic_peaks()
                    # if abs(precursor_mz - spec.prec_mz) < 0.05:
                    if abs(precursor_mz - spec.prec_mz) < 0.05:
                        query_spec_list.append(spec)
                    # except:
                    # print('mzML converted error')
    print("read {} mzml file, {} spectra complete".format(len(file_list), len(query_spec_list)))

    # READ query block and SEARCH against bmdms
    list_bmdms_inchikey = repo_spectrum.get_inchikey_list_metfrag_valid(
        instrument_type=instrument_type, precursor_type=ion_mode
    )
    print("bmdms inchikey #: {}\n".format(len(list_bmdms_inchikey)))
    print("bmdms_spectrum_dict_preparing... {}\n".format(datetime.datetime.now()))
    bmdms_spectrum_dict = repo_spectrum.get_spectrum_dict_by_inchikey_choose_metfrag_valid(
        list_bmdms_inchikey,
        instrument_type=instrument_type,
        precursor_type=ion_mode,
        include_len_0=False,
        sample_num=20,
    )

    repo_spectrum.close()
    
    list_bmdms_inchikey = list(bmdms_spectrum_dict.keys())
    print("bmdms_spectrum_dict_prepared {}\n".format(datetime.datetime.now()))
    print("bmdms inchikey #: {}\n".format(len(list_bmdms_inchikey)))

    print("Noise mod = {}\n".format(NOISE_MOD))

    # print(
    #     "\ttarget, bmdms spec length {}\n".format(
    #         [
    #             len(bmdms_spectrum_dict[bmdms_inchikey])
    #             for bmdms_inchikey in list_bmdms_inchikey
    #         ]
    #     )
    # )

    template, score_compare, score_count = [], [], []

    for target_inchikey in list_bmdms_inchikey:

        target_spec_list = bmdms_spectrum_dict[target_inchikey]

        count = 0
        for target_spec in target_spec_list:
            for query_spec in query_spec_list:
                target_peaks = target_spec.peaks_cut_noise(mod=NOISE_MOD).tolist()
                query_peaks = query_spec.peaks
                if len(target_peaks) > 4 and len(query_peaks) > 4:
                    template.append((target_peaks, query_peaks))
                    count += 1

        score_count.append([target_inchikey, count])

    print("\t\t\tPairs to calculate : {}\n".format(len(template)))

    for x1 in tool_list.divide_list(template, MULTI_PROC_UNIT):
        pool = multiprocessing.Pool(processes=NUM_PROCESSOR)
        score_compare = score_compare + pool.starmap(
            tool_calculate.SimilarityCalculator.calculate_cosine_correlation,
            x1,
        )
        pool.close()
        pool.join()

    inchikey_temp, score_temp = tool_calculate.divide_count(
        score_count=score_count, score_compare=score_compare
    )
    
    rank_index_array = tool_list.get_index_value_top_k(np.array(score_temp), 20)
    for n, rank_index in enumerate(rank_index_array):
        # if score_temp[rank_index] > 0.5:
        print("{}th\tinchikey = {}\tscore={}".format(n+1, inchikey_temp[rank_index], score_temp[rank_index]))

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print(main.__doc__)
        exit()
    MZML_FOLDER = sys.argv[1]
    INSTRUMENT_TYPE = sys.argv[2]
    if INSTRUMENT_TYPE in ('Orbitrap', 'O'):
        INSTRUMENT_TYPE = 'Orbitrap'
    elif INSTRUMENT_TYPE in ('QTOF', 'Q-TOF', 'Q'):
        INSTRUMENT_TYPE = 'QTOF'
    else:
        print('Instrument type error')
        sys.exit()
    ION_MODE = sys.argv[3]
    if ION_MODE in ('[M+H]+', 'H'):
        ION_MODE = '[M+H]+'
    elif ION_MODE in ('[M+Na]+', 'Na'):
        ION_MODE = '[M+Na]+'
    else:
        print('ion mode(precursor ion) type error')
        sys.exit()
    PRECURSOR_MZ = float(sys.argv[4])
    main(query_path_dir=MZML_FOLDER, instrument_type=INSTRUMENT_TYPE, ion_mode=ION_MODE, precursor_mz=PRECURSOR_MZ)
