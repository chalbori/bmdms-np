from decimal import Decimal
import copy
import multiprocessing
import numpy as np
from scipy import spatial
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs
from tools.tool_peaks import merge_peaklist, delete_zero_intensity


class SimilarityCalculator:
    @staticmethod
    def multi_cal_cosine_correlation(num_processor, exp_spec_list, public_spec_list):

        pool = multiprocessing.Pool(processes=4)
        template_list = []
        for exp_spec in exp_spec_list:
            for pub_spec in public_spec_list:
                template_list.append((exp_spec.peaks.tolist(), pub_spec.peaklist))
        score_list = pool.starmap(
            SimilarityCalculator.calculate_cosine_correlation, template_list
        )
        pool.close()
        pool.join()
        return score_list

    @staticmethod
    def calculate_cosine_correlation(query_peaklist, target_peaklist, mod="align"):

        if mod == "align":
            query_input_peaklist, target_input_peaklist = SimilarityCalculator.align_for_similarity(
                query_peaklist, target_peaklist
            )
        # elif mod == 'bin':
        #     query_str_peaklist = [[str(x[0]), str(x[1])] for x in query_peaklist]
        #     target_str_peaklist = [[str(x[0]), str(x[1])] for x in target_peaklist]
        #     to_num_decimal_place = max(_num_decimal_places(query_str_peaklist), _num_decimal_places(target_str_peaklist))
        #     query_input_peaklist = Peaklist.bin_peaklist(to_num_decimal_place, query_peaklist, mode='max')
        #     target_input_peaklist = Peaklist.bin_peaklist(to_num_decimal_place, target_peaklist, mode='max')

        # elif mod == 'bin_zero':
        #     query_str_peaklist = [[str(x[0]), str(x[1])] for x in query_peaklist]
        #     target_str_peaklist = [[str(x[0]), str(x[1])] for x in target_peaklist]
        #     to_num_decimal_place = max(_num_decimal_places(query_str_peaklist),
        #                                _num_decimal_places(target_str_peaklist))
        #     query_input_peaklist = Peaklist.bin_peaklist(0, query_peaklist, mode='max')
        #     target_input_peaklist = Peaklist.bin_peaklist(0, target_peaklist, mode='max')

        elif mod == "":
            query_input_peaklist = query_peaklist
            target_input_peaklist = target_peaklist

        q_peaklist = np.array(
            query_input_peaklist, dtype="float64"
        )  # prep: np array initialization + weighting
        t_peaklist = np.array(target_input_peaklist, dtype="float64")

        q_peaklist_weight = SimilarityCalculator.weight(q_peaklist)
        t_peaklist_weight = SimilarityCalculator.weight(t_peaklist)

        comb1 = SimilarityCalculator.combine(
            pl1=q_peaklist_weight, pl2=t_peaklist_weight
        )
        comb2 = SimilarityCalculator.combine(
            pl1=t_peaklist_weight, pl2=q_peaklist_weight
        )
        return 1 - spatial.distance.cosine(comb1[:, 1], comb2[:, 1])

    @staticmethod
    def align_for_similarity(pl1, pl2):
        # from ms2chem.data.peaklist import align, cluster
        list_np_peaklist = [np.array(pl, dtype="float64") for pl in [pl1, pl2]]
        merged_peaklist = merge_peaklist(list_np_peaklist)
        merged_peaklist = delete_zero_intensity(merged_peaklist)
        merged_aligned_peaklist = align(merged_peaklist.tolist(), 5, -2)
        merged_cluster = cluster(merged_peaklist.tolist(), 5, -2)

        dict_pl1 = {peak[0]: peak[1] for peak in pl1}
        dict_pl2 = {peak[0]: peak[1] for peak in pl2}
        dict_pl1_keys = dict_pl1.keys()
        dict_pl2_keys = dict_pl2.keys()

        cluster_peak1 = copy.deepcopy(merged_cluster)
        cluster_peak2 = copy.deepcopy(merged_cluster)

        for bin in cluster_peak1:
            for peak in bin:
                if peak[0] in dict_pl1_keys:
                    peak[1] = dict_pl1.get(peak[0])
                else:
                    bin.remove(peak)

        for bin in cluster_peak2:
            for peak in bin:
                if peak[0] in dict_pl2_keys:
                    peak[1] = dict_pl2.get(peak[0])
                else:
                    bin.remove(peak)

        result_pl1 = copy.deepcopy(merged_aligned_peaklist)
        result_pl2 = copy.deepcopy(merged_aligned_peaklist)

        for i in range(0, len(result_pl1)):
            if cluster_peak1[i] != []:
                result_pl1[i][1] = max(intensity[1] for intensity in cluster_peak1[i])
            else:
                result_pl1[i][1] = 0

        for i in range(0, len(result_pl2)):
            if cluster_peak2[i] != []:
                result_pl2[i][1] = max(intensity[1] for intensity in cluster_peak2[i])
            else:
                result_pl2[i][1] = 0

        return result_pl1, result_pl2

    @staticmethod
    def combine(pl1, pl2):
        copied_pl2 = pl2.copy()
        copied_pl2[:, 1] = 0
        temp = np.concatenate((pl1, copied_pl2), axis=0)
        _, i = np.unique(temp[:, 0], return_index=True)
        combine_unique = temp[np.sort(i)]
        return combine_unique[combine_unique[:, 0].argsort()]

    @staticmethod
    def weight(peaklist, weight_mtoz=2, weight_intensity=0.5):
        return np.array(
            [
                peaklist[:, 0],
                np.power(peaklist[:, 0], weight_mtoz)
                * np.power(peaklist[:, 1], weight_intensity),
            ]
        ).T


def compare_structure(smiles1, smiles2, fp_type="Morgan", sim_type="Dice"):
    """
    Task: Compare structual similarity of two compound based on fingerprints.
    Parameters:
        smiles1: str, smiles of the compound 1
        smiles2: str, smiles of the compound 2
        fp_type: str, type of fingerprints
        sim_type: str, method for calculating similarity
    """
    if fp_type == "Morgan":
        getfp = lambda smi: AllChem.GetMorganFingerprint(
            Chem.MolFromSmiles(smi), 2, useFeatures=False
        )
    elif fp_type == "MorganWithFeature":
        getfp = lambda smi: AllChem.GetMorganFingerprint(
            Chem.MolFromSmiles(smi), 2, useFeatures=True
        )
    elif fp_type == "MACCS":
        getfp = lambda smi: Chem.MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smi))
    elif fp_type == "Topological":
        getfp = lambda smi: FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smi))
    elif fp_type == "AtomPairs":
        getfp = lambda smi: Pairs.GetAtomPairFingerprint(Chem.MolFromSmiles(smi))

    try:
        fp1 = getfp(smiles1)
        fp2 = getfp(smiles2)
        if sim_type == "Dice":
            sim_fp = DataStructs.DiceSimilarity(fp1, fp2)
        elif sim_type == "Tanimoto":
            sim_fp = DataStructs.TanimotoSimilarity(fp1, fp2)
        elif sim_type == "Cosine":
            sim_fp = DataStructs.CosineSimilarity(fp1, fp2)
        elif sim_type == "Sokal":
            sim_fp = DataStructs.SokalSimilarity(fp1, fp2)
        elif sim_type == "Russel":
            sim_fp = DataStructs.RusselSimilarity(fp1, fp2)

    except Exception as e:
        sim_fp = -1
    return sim_fp


# m/z 기준으로 더해진 peaklist, error range : +/- error_value * 10 ^ (-num_decimal_places)
def align(merged_peaklist, error_value, num_decimal_places):
    clustered_peaklist = cluster(merged_peaklist, error_value, num_decimal_places)

    # if [[]] in binned_peaklist:
    result = []
    for x in clustered_peaklist:
        result.append(get_weighted_centroid(x))
    return copy.deepcopy(result)


def cluster(merged_peaklist, error_value, num_decimal_places):
    binned_peaklist = read_bin(merged_peaklist, error_value, num_decimal_places)
    if len(binned_peaklist) < 2:
        return [copy.deepcopy(merged_peaklist)]
    else:
        for i in range(0, len(binned_peaklist) - 1):
            # print(binned_peaklist[i], binned_peaklist[i+1])
            if abs(binned_peaklist[i][-1][0] - binned_peaklist[i + 1][0][0]) < (
                2 * error_value * pow(10, num_decimal_places)
            ):
                result = hierarchy_cluster(
                    binned_peaklist[i],
                    binned_peaklist[i + 1],
                    error_value,
                    num_decimal_places,
                )
                # print(i, result)
                # print(type(result[0]), type(result[1]))

                binned_peaklist[i] = result[0]
                binned_peaklist[i + 1] = result[1]
            else:
                pass
        result = []
        for x in binned_peaklist:
            if x != [[]]:
                result.append(x)
                # binned_peaklist.remove([[]])
    # return binned_peaklist
    return result


def get_weighted_centroid(peaklist):
    sum_m_i = 0.0
    sum_i = 0.0
    # print('input_calculate centroid:',peaklist)
    for i in range(0, len(peaklist)):
        sum_m_i += peaklist[i][0] * peaklist[i][1]
        sum_i += peaklist[i][1]
    return [sum_m_i / sum_i, max([x[1] for x in peaklist])]


def _num_decimal_places(formatted_peaklist):
    # type check. formatted_peaklist must be a list of list(str)
    if not issubclass(type(formatted_peaklist[0][0]), str):
        raise TypeError(
            "{!r} is not a subclass of {!r}".format(type(formatted_peaklist), str)
        )

    min_dp = 5
    dp_cnt_generator = (
        Decimal(peak[0]).as_tuple().exponent for peak in formatted_peaklist
    )
    dp_cnt_dict = dict()
    all_dp_cnt = 0
    for dp in dp_cnt_generator:
        all_dp_cnt += 1
        each_dp_cnt = dp_cnt_dict.get(dp, None)
        if each_dp_cnt:
            dp_cnt_dict[each_dp_cnt] += 1
        else:
            dp_cnt_dict.update({dp: 0})

    dp_percent_dict = {dp: dp_cnt_dict.get(dp) / all_dp_cnt for dp in dp_cnt_dict}

    for dp in dp_percent_dict.keys():
        if dp_percent_dict.get(dp) >= 0.8:
            return dp
        elif dp < min_dp:
            min_dp = dp
    return min_dp


def read_bin(peaklist, error_value, num_decimal_places):
    if not len(peaklist) > 1:
        return [peaklist.copy()]

    temp_peaklist = peaklist.copy()
    new_bin = []

    new_bin_temp = [temp_peaklist[0]]
    temp_peaklist.remove(temp_peaklist[0])

    for peak in temp_peaklist:
        if peak[0] - new_bin_temp[0][0] < 2 * error_value * pow(10, num_decimal_places):
            new_bin_temp.append(peak)
            # temp_peaklist.remove(peak)
        else:
            new_bin.append(new_bin_temp)
            new_bin_temp = [peak]
            # temp_peaklist.remove(peak)
    new_bin.append(new_bin_temp)
    return new_bin


def hierarchy_cluster(binned_a, binned_b, error_value, num_decimal_places):
    from scipy.cluster.hierarchy import linkage
    import collections

    bin_size = 2 * error_value * pow(10, num_decimal_places)
    # 논문 상의 2E = bin_size
    # if abs(binned_a[-1][0] - binned_b[0][0]) < 2 * pow(10, num_decimal_places):
    temp_list = copy.deepcopy(binned_a + binned_b)
    index_list = []
    for i in range(0, len(temp_list)):
        index_list.append([int(i)])
    # templist = [[x[0] / bin_size, x[1]] for x in temp_list]
    z = linkage(temp_list, "complete")

    cluster_dict = dict()
    for i in range(0, len(z)):
        new_cluster_id = len(temp_list) + i  # templist?
        old_cluster_id_0 = int(z[i, 0])
        old_cluster_id_1 = int(z[i, 1])
        combined_ids = list()
        if old_cluster_id_0 in cluster_dict:
            combined_ids += cluster_dict[old_cluster_id_0]
            # del cluster_dict[old_cluster_id_0]
        else:
            combined_ids += [old_cluster_id_0]
        if old_cluster_id_1 in cluster_dict:
            combined_ids += cluster_dict[old_cluster_id_1]
            # del cluster_dict[old_cluster_id_1]
        else:
            combined_ids += [old_cluster_id_1]
        cluster_dict[new_cluster_id] = combined_ids
    cluster_dict_ordered = collections.OrderedDict(sorted(cluster_dict.items()))

    if int(z[-1][1]) >= len(temp_list):
        cluster_last_index_list = cluster_dict_ordered[int(z[-1][1])]
        cluster_last_list = [temp_list[x] for x in cluster_last_index_list]
    else:
        cluster_last_list = [temp_list[int(z[-1][1])]]
    sorted_cll = sorted(cluster_last_list, key=lambda l: l[0])

    if int(z[-1][0]) >= len(temp_list):
        cluster_second_last_index_list = cluster_dict_ordered[int(z[-1][0])]
        cluster_second_last_list = [
            temp_list[x] for x in cluster_second_last_index_list
        ]
    else:
        cluster_second_last_list = [temp_list[int(z[-1][0])]]
    sorted_csll = sorted(cluster_second_last_list, key=lambda l: l[0])

    # if z[-1][2] >= 1: # if distance between clusters >= 1 return each cluster, each cluster will be calculated with centroid
    #     if sorted_csll[-1][0] < sorted_cll[0][0]:
    #         return sorted_csll, sorted_cll
    #     else:
    #         return sorted_cll, sorted_csll
    # else: # if distance between cluster < 1 return combined cluster. The cluster will be calculated centroid
    #     return [[]], sorted(sorted_csll + sorted_cll, key=lambda l: l[0])

    if (
        abs(
            get_weighted_centroid(sorted_csll)[0] - get_weighted_centroid(sorted_cll)[0]
        )
        < bin_size
    ):
        if sorted_csll[-1][0] < sorted_cll[0][0]:
            return [[]], sorted_csll + sorted_cll
        else:
            return [[]], sorted_cll + sorted_csll

    else:
        if sorted_csll[-1][0] < sorted_cll[0][0]:
            return sorted_csll, sorted_cll
        else:
            return sorted_cll, sorted_csll


def divide_count(score_count, score_compare):
    y_list = []
    max_score_list = []
    for score in score_count:
        score_list = score_compare[: score[1]]
        if len(score_list) > 0:
            y_list.append(score[0])
            max_score_list.append(max(score_list))
            score_compare = score_compare[score[1] :]
    if len(score_compare) != 0:
        print("Error in score_list")
    return y_list, max_score_list
