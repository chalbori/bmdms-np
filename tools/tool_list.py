import numpy as np


def divide_list(l, n):
    """
    yield successive n-sized
    :param l:
    :param n:
    :return:
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_overlapped_set(query_set, target_set, limit_count):
    result_set = set()
    for query in query_set:
        if query in target_set:
            result_set.add(query)
        if len(result_set) == limit_count:
            break
    return result_set


def get_index_value_top_k(nparray, k):
    if not any(nparray>0): # avoid nparray=([0.0, 0.0, ..., 0.0])
        print("All score is 0.0")
        # return np.zeros(1)
        return [0]
    if np.unique(nparray[nparray.argsort()[-k:][::-1]]).size == 1:
        return np.where(nparray == nparray.max())[0]
    else:
        return nparray.argsort()[-k:][::-1]
