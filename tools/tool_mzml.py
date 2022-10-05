def extract_prec_mz(mzml):
    for precursor_dict in mzml.selected_precursors:
        # print('precursor m/z : ', precursor_dict['mz'])
        return precursor_dict['mz']

from tools import tool_file
from pyopenms import MSExperiment, MzMLFile


def read_mzml(input_folder: str):
    ms_exp_list = []
    file_list = tool_file.list_dir_except_hidden(path=input_folder, suffix=".mzML")
    for mzml_file in file_list:
        print("file read: {}".format(mzml_file))
        exp = MSExperiment()
        MzMLFile().load(mzml_file, exp)
        ms_exp_list.append(exp)
    return ms_exp_list
