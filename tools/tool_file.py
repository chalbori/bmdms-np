import os
import csv
from tools import tool_chemical


def list_dir_except_hidden(path, suffix="", read_sub=False):
    file_list = []
    if read_sub:
        for path, _, files in os.walk(path):
            for name in files:
                if not name.startswith('.'):
                    if name.endswith(suffix):
                        file_list.append(os.path.join(path, name))
    else:
        for file_name in os.listdir(path):
            if not file_name.startswith('.'):
                if file_name.endswith(suffix):
                    file_list.append(os.path.join(path, file_name))
    file_list.sort()
    return file_list


def get_year_batch_group_ce_info(file):
    filename = os.path.splitext(file)[0]
    y_b_g = filename.split('_')
    return int(y_b_g[1]), int(y_b_g[3]), int(y_b_g[5]), int(y_b_g[7])


def read_tsv(file):
    """
    read and remove duplicate str
    :param file:
    :return: list
    """
    f = open(file=file, encoding='utf-8')
    rdr = csv.reader(f, delimiter='\t')
    str_list = set()
    for line in rdr:
        str_list.add(line[0])
    f.close()
    return list(str_list)


def read_smi(file):
    """
    read and remove duplicate str
    :param file:
    :return: list
    """
    f = open(file=file, encoding='utf-8')
    rdr = csv.reader(f, delimiter='\t')
    str_list = set()
    for line in rdr:
        # str_list.add(cancel_isomeric_smi(line[0]))
        str_list.add(line[0])
    f.close()
    return list(str_list)


def cancel_isomeric_smi(smi_string):
    import data.chemical as chemical
    from rdkit.Chem import AllChem
    smiles = ""
    try:
        mol = chemical.read_string("smi", smi_string)
        smiles = AllChem.MolToSmiles(mol, isomericSmiles=False)
    except:
        print(smi_string, "Error in conversion")
    return smiles


def convert_smi_to_inchikey(smi_string):
    import data.chemical as chemical
    from rdkit.Chem import AllChem
    smiles = ""
    try:
        mol = chemical.read_string("smi", smi_string)
        inchikey = AllChem.MolToInchiKey(mol)
    except:
        print(smi_string, "Error in conversion")
    return inchikey


# def mkdir_folder(directory):
#     import os
#     if not os.path.exists(directory):
#         os.makedirs(directory)
def mkdir_folder(directory):
    import os, errno
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def dump_pickle_file(file, contents):
    import gzip, pickle
    with open(file, 'wb') as f:
        pickle.dump(contents, f)


def read_pickle_file_gz(file):
    import gzip, pickle
    with open(file, 'rb') as f:
        result = pickle.load(f)
    return result


def dump_json_file(file, contents):
    import json
    with open(file, 'w', encoding='utf-8') as f:
        json.dump(contents, f)


def read_json_file(file):
    import json
    with open(file, 'r', encoding='utf-8') as f:
        result = json.load(f)
    return result


def read_top_k(file):
    """

    :param file:
    :return:
    """
    f = open(file=file, encoding='utf-8')
    rdr = csv.DictReader(f, delimiter=',')
    top_k, list_true, list_false, list_tp = [], [], [], []
    for row in rdr:
        top_k.append(int(row['top']))
        list_true.append(int(row['TRUE']))
        list_false.append(int(row['FALSE']))
        list_tp.append(int(row['TRUE']) / (int(row['TRUE']) + int(row['FALSE'])) * 100)
    f.close()
    return top_k, list_true, list_false, list_tp


def convert_spec_to_msp(chemical, spec, peaklist):
    result = ""
    result += 'NAME: ' + chemical.name + '\n'
    result += 'PRECURSORMZ: ' + str(spec.prec_mz) + '\n'
    result += 'PRECURSORTYPE: ' + spec.prec_type + '\n'
    result += 'FORMULA: ' + chemical.molecular_formula + '\n'
    result += 'Ontology: ' + chemical.name + '\n'
    result += 'INCHIKEY: ' + spec.inchikey + '\n'
    result += 'SMILES: ' + spec.smiles + '\n'
    result += 'RETENTIONTIME: ' + '' + '\n' 
    result += 'CCS: ' + '' + '\n'
    result += 'IONMODE: ' + 'Positive' + '\n'
    result += 'INSTRUMENTTYPE: ' + spec.instrument_type + '\n'
    result += 'INSTRUMENT: ' + '' + '\n'
    result += 'COLLISIONENERGY: ' + str(spec.collision_energy) + '\n'
    result += 'Comment: ' + 'spec_id=' + str(spec.spec_id) + '; ' + 'origin=BMDMS-NP' + '\n'
    result += 'Num Peaks: ' + str(len(peaklist)) + '\n'
    for peak in peaklist:
        result += str(peak[0]) + '\t' + str(peak[1]) + '\n'

    """
    ref: http://prime.psc.riken.jp/compms/msdial/download/msp/MSMS-Pos-MassBankEU.msp
    NAME: CLC_301.1468_14.3
    PRECURSORMZ: 301.14660285217
    PRECURSORTYPE: [M+H]+
    FORMULA: C18H21ClN2
    Ontology: Diphenylmethanes
    INCHIKEY: WFNAKBGANONZEQ-UHFFFAOYNA-N
    SMILES: CN1CCN(CC1)C(C2=CC=CC=C2)C3=CC=C(C=C3)Cl
    RETENTIONTIME: 
    CCS: 178.3135852
    IONMODE: Positive
    INSTRUMENTTYPE: LC-ESI-QFT
    INSTRUMENT: Q Exactive Orbitrap Thermo Scientific
    COLLISIONENERGY: 15, 30, 45, 60, 70 or 90 (nominal)
    Comment: DB#=ET010001; origin=MassBank-EU
    Num Peaks: 3
    165.0694	20
    166.0777	30
    201.0466	1000
    """
    return result
