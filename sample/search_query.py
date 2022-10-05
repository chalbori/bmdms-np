import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
print(os.path.abspath(os.path.dirname(__file__)))
import io
import data.spectrum as spectrum
from database import info_chem_mysql, info_spectrum_mysql
import pymzml
import numpy as np
import sqlite3
import tools.tool_file as tool_file
import tools.tool_mzml as tool_mzml
import tools.tool_list as tool_list
from datetime import datetime
import multiprocessing
import tools.tool_calculate as tool_calculate

# Asiaticoside      c20150100009s
# Madecassoside     c20160300306
# Asiaticoside B    c20160300806
# Asiatic acid      c20150300081
# Madecassic acid   c20150300186
# Magnolol          c20150100069
# Honokiol          c20150100059
target_h_dict = {'c20150100009': 959.5216,
                 'c20160300306': 975.5165,
                 'c20160300806': 975.5165,
                 'c20150300081': 489.3581,
                 'c20150300186': 505.3530,
                 'c20150100069': 267.1385,
                 'c20150100059': 267.1385}
target_na_dict = {'c20150100009': 981.5034,
                  'c20160300306': 997.4984,
                  'c20160300806': 997.4984,
                  'c20150300081': 511.3400,
                  'c20150300186': 527.3349,
                  'c20150100069': 289.1204,
                  'c20150100059': 289.1204}

if __name__ == '__main__':
    num_processor = 4
    OUT_FOLDER = 'out/fig/query_search_orbitrap/'
    tool_file.mkdir_folder(OUT_FOLDER)
    file_exp_smiles_list = 'resources/expdb_smiles_list.pickle.gz'
    smi_set_count = 1000
    multi_proc_unit = 2000
    top_k_num = 20
    # noise_mod = ['original', 'mean', 'median', 'dynamic']
    noise_mod = ['dynamic']
    ion_mode = '[M+H]+'
    instrument_type_expdb = 'Orbitrap'
    log_file = OUT_FOLDER + 'log_{}.txt'.format(datetime.now())
    log = open(log_file, 'w')

    db_exp_uri = 'sqlite:///resources/ms2db_orbitrap_qtof_all.sqlite'
    repo_chemical = info_chem_mysql.CompoundRepository(db_uri=db_exp_uri)
    repo_spectrum = info_spectrum_mysql.SpectrumRepository(db_uri=db_exp_uri)
    query_path_dir = 'resources/20180903_Biospectrume_Orbitrap/20180903_Extract_mzML/Centellaasiatica_1/'
    file_list = tool_file.list_dir_except_hidden(path=query_path_dir)
    file_list.sort()

    spec_list = []
    for mzml_file in file_list:
        run = pymzml.run.Reader(mzml_file)

        # prec_list = []
        for n, mzml in enumerate(run):
            if mzml.ms_level == 2 and len(mzml.peaks('centroided')) > 4:
                try:
                    spec = spectrum.Spectrum()
                    spec.chem_id = None
                    spec.c_id= None
                    spec.prec_mz = tool_mzml.extract_prec_mz(mzml=mzml)
                    # prec_list.append(spec.prec_mz)
                    spec.instrument_type = 'Orbitrap'
                    # spec.collision_energy = int(collision_energy)
                    spec.peaks = mzml.peaks('centroided')
                    spec.peaks = spec.peaks[spec.peaks[:, 1] != 0]
                    spec.noise_mean = mzml.estimated_noise_level('mean')
                    spec.noise_median = mzml.estimated_noise_level('median')
                    spec.noise_mad = mzml.estimated_noise_level('mad')
                    spec.noise_dynamic = spec.estimate_noise_dynamic()
                    spec_list.append(spec)
                except:
                    print('mzML converted error')
    # print(len(spec_list))

    query_dict = {}
    for moleucle in target_h_dict.keys():
        temp_list = []
        # print(target_h_dict[moleucle])
        for spec in spec_list:
            # print('\t',spec.prec_mz)
            if abs(target_h_dict[moleucle] - spec.prec_mz) < target_h_dict[moleucle] * 100 / 1000000:
                # print(target_h_dict[moleucle], spec.prec_mz)
                temp_list.append(spec)
        if len(temp_list):
            query_dict[moleucle] = temp_list

    for key in query_dict.keys():
        log.write('c_id: {}\tlen(query_list){}\n'.format(key, len(query_dict[key])))

    query_keys = list(query_dict.keys())

    # To READ and CALCULATE
    expdb_cid_list = repo_spectrum.get_cid_list(instrument_type=instrument_type_expdb, precursor_type=ion_mode)
    log.write('exp ids #: {}\n'.format(len(expdb_cid_list)))
    log.flush()
    log.write('exp_spectrum_dict_preparing... {}\n'.format(datetime.now()))
    exp_spec_dict = repo_spectrum.get_spectrum_dict_by_cid_choose(expdb_cid_list,
                                                                  instrument_type=instrument_type_expdb,
                                                                  precursor_type=ion_mode, include_len_0=False,
                                                                  sample_num=15)

    expdb_cid_list = list(exp_spec_dict.keys())
    exp_smi_precmz_dict = {}
    for cid in expdb_cid_list:
        expdb_prec_mz_list = [spec.prec_mz for spec in exp_spec_dict[cid]]
        if len(expdb_prec_mz_list):
            exp_smi_precmz_dict[cid] = expdb_prec_mz_list
    log.flush()

    # noise = 'dynamic'
    for noise in noise_mod:
        log.write('Noise mod = {}\n'.format(noise))
        y, score, top_k = [], [], []  #
        # y_na, score_na, top_k_na = [], [], []  #
        # y_both, score_both, top_k_both = [], [], []

        log.write('\ttarget, exp spec length {}\n'.format(
            [len(exp_spec_dict[target_smiles]) for target_smiles in expdb_cid_list]))

        total_idx = 0
        for n in range(len(query_keys)):
            query_id = query_keys[n]
            log.write(
                '{} / {} {} calculation of similarity processed {}. {}\n'.format(n + 1, len(query_keys),
                                                                                 query_id,
                                                                                 datetime.now(), noise))
            top_k_h_temp, template_h, score_compare, score_count = [], [], [], []

            query_spec_list = query_dict[query_id]
            log.write('\tquery, public spec length {}\n'.format(len(query_spec_list)))

            query_spec_prec_mz_list = [spec.prec_mz for spec in query_spec_list]

            min_prec_mz = float(min(query_spec_prec_mz_list))
            max_prec_mz = float(max(query_spec_prec_mz_list))

            target_cid_list = []
            for c_id, prec_mz_list in exp_smi_precmz_dict.items():
                if float(min_prec_mz - 0.1) < min(prec_mz_list):
                    if float(max_prec_mz + 0.1) > max(prec_mz_list):
                        target_cid_list.append(c_id)

            for cid in target_cid_list:
                target_spec_list = exp_spec_dict[cid]

                count_h = 0
                for exp_spec in target_spec_list:
                    for query_spec in query_spec_list:
                        exp_peaks = exp_spec.peaks_cut_noise(mod=noise).tolist()
                        pub_peaks = query_spec.peaks
                        if len(exp_peaks) > 4 and len(pub_peaks) > 4:
                            template_h.append((exp_peaks, pub_peaks))
                            count_h += 1

                score_count.append([cid, count_h])

            log.write('\t\t\tPairs to calculate : {}\n'.format(len(template_h)))
            log.flush()

            for x1 in tool_list.divide_list(template_h, multi_proc_unit):
                pool = multiprocessing.Pool(processes=num_processor)
                score_compare = score_compare + pool.starmap(
                    tool_calculate.SimilarityCalculator.calculate_cosine_correlation, x1)
                pool.close()
                pool.join()

            y_temp, score_temp = tool_calculate.divide_count(score_count=score_count,
                                                             score_compare=score_compare)

            pickle_y_temp = OUT_FOLDER + 'y_{}_{}th.pickle.gz'.format(noise, total_idx)
            pickle_score_temp = OUT_FOLDER + 'score_{}_{}th.pickle.gz'.format(noise, total_idx)
            tool_file.dump_pickle_file(pickle_y_temp, y_temp)
            tool_file.dump_pickle_file(pickle_score_temp, score_temp)
            total_idx += 1


def extract_prec_mz(mzml):
    for precursor_dict in mzml.selected_precursors:
        # print('precursor m/z : ', precursor_dict['mz'])
        return precursor_dict['mz']


def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)
