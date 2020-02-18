import numpy as np


def _all_same(items):
    return all(x == items[0] for x in items)


def merge_peaklist(list_np_peaklist):
    temp = np.zeros((1, 2))
    for np_peaklist in list_np_peaklist:
        temp = sum(temp, np_peaklist)

    temp = normalize(temp)

    return temp


def sum(np_pl1, np_pl2):
    combined_1 = combine(np_pl1, np_pl2)
    combined_2 = combine(np_pl2, np_pl1)
    temp = np.append([combined_1[:, 0]],
                     [np.maximum(combined_1[:, 1], combined_2[:, 1])], axis=0).T
    return temp


def normalize(np_peaklist):
    max_intensity = np.max(np_peaklist[:, 1])
    np_peaklist[:, 1] /= max_intensity
    return np_peaklist


def combine(pl1, pl2):
    copied_pl2 = pl2.copy()
    copied_pl2[:, 1] = 0
    temp = np.concatenate((pl1, copied_pl2), axis=0)
    _, i = np.unique(temp[:, 0], return_index=True)
    combine_unique = temp[np.sort(i)]
    return combine_unique[combine_unique[:, 0].argsort()]


def delete_zero_intensity(np_peaklist):
    return np_peaklist[np.all(np_peaklist != 0, axis=1)]
