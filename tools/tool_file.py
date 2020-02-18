import os
import csv


def list_dir_except_hidden(path, read_sub=False):
    file_list = []
    if read_sub:
        for path, _, files in os.walk(path):
            for name in files:
                if not name.startswith('.'):
                    file_list.append(os.path.join(path, name))
    else:
        for file_name in os.listdir(path):
            if not file_name.startswith('.'):
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
