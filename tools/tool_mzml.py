def extract_prec_mz(mzml):
    for precursor_dict in mzml.selected_precursors:
        # print('precursor m/z : ', precursor_dict['mz'])
        return precursor_dict['mz']