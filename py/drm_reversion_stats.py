import re

import numpy as np
import pandas as pd
import surpyval
from pastml import datetime2numeric
from pastml.tree import read_forest, DATE, annotate_dates

MIN_CENSORED_DATA_POINT_NUMBER = 5

def split_drm(drm):
    """
    Splits a given DRM into (prefix, to_letter)
    :param drm: input DRM
    :return: tuple (prefix, to_letter)
    """
    res = re.findall(r'((RT|PR|IN)[:_][A-Z][0-9]+)([A-Z]+)', drm)
    return (res[0][0], res[0][-1]) if res else (None, None)


def survival_estimate(loss_durations, censorship):
    if not len(loss_durations) or sum(((censorship == -1) | (censorship == 2)).astype(int)) < MIN_CENSORED_DATA_POINT_NUMBER \
            or sum(((censorship == 1) | (censorship == 2)).astype(int)) < MIN_CENSORED_DATA_POINT_NUMBER:
        return None

    loss_est = None

    min_ld, max_ld, sum_ld = np.inf, -np.inf, 0
    for ld in loss_durations:
        if isinstance(ld, list):
            ld = (ld[0] + ld[1]) / 2
        if ld < min_ld:
            min_ld = ld
        if ld > max_ld:
            max_ld = ld
        sum_ld += ld

    for i in range(10):
        if i == 0:
            x = sum_ld / len(loss_durations)
        else:
            x = np.random.uniform(min_ld, max_ld)
        model = surpyval.Weibull.fit(x=loss_durations, c=censorship, fixed={'beta': 1}, init=[x, 1])
        if model.res.success:
            return model.params[0]
        else:
            print(loss_durations)
            print(censorship)
            exit()
    return loss_est


def get_arv_class(arv_df, prefix, letter):
    sub_arv_df = arv_df[(arv_df['prefix'] == prefix) & arv_df['letters'].apply(lambda _: letter in _)]
    return next(iter(sub_arv_df['drug class'].unique()))


def estimate_loss_time_and_ci(loss_durations, censorship, repetitions=1000):
    loss_est_combo = survival_estimate(loss_durations.tolist(), censorship)
    ci_min, ci_max = loss_est_combo, loss_est_combo

    if loss_est_combo is not None:
        indices = np.array(range(len(loss_durations)))
        bootstrapped_estimates = []
        for i in range(repetitions):
            indices_bootstrapped = np.hstack([
                np.random.choice(indices[(censorship == -1) | (censorship == 2)], MIN_CENSORED_DATA_POINT_NUMBER, replace=True),
                np.random.choice(indices[(censorship == 1) | (censorship == 2)], MIN_CENSORED_DATA_POINT_NUMBER, replace=True),
                np.random.choice(indices, len(indices) - 2 * MIN_CENSORED_DATA_POINT_NUMBER, replace=True)])

            bootstrapped_estimates.append(survival_estimate(loss_durations[indices_bootstrapped].tolist(),
                                                            censorship[indices_bootstrapped]))
        bootstrapped_estimates = np.array(bootstrapped_estimates)
        ci_min = np.quantile(bootstrapped_estimates, 0.025)
        ci_max = np.quantile(bootstrapped_estimates, 0.975)

    return loss_est_combo, ci_min, ci_max


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--drms', type=str)
    parser.add_argument('--poly', type=str)
    parser.add_argument('--drm_data', type=str)
    parser.add_argument('--treatment_data', type=str)
    parser.add_argument('--Castro_data', type=str)
    parser.add_argument('--arv_data', type=str)
    parser.add_argument('--mp_pattern', type=str)
    parser.add_argument('--nwk', type=str)
    parser.add_argument('--output_tab', type=str)
    parser.add_argument('--subtype', type=str)
    params = parser.parse_args()

    poly_df = pd.read_csv(params.poly, sep='\t', index_col=0)

    with open(params.drms, 'r') as f:
        drms_by_pos = sorted(f.read().strip().split(' '))

    drms = []
    for _ in drms_by_pos:
        prefix, letters = split_drm(_)
        drms.extend('{}{}'.format(prefix, l) for l in letters)
    drm_df = pd.read_csv(params.drm_data, sep='\t', index_col=0)[drms]
    drm_df.index = drm_df.index.map(str)

    tree = read_forest(params.nwk)[0]
    annotate_dates([tree])
    tree_ids = {_.name for _ in tree.traverse()}

    treatment_df = pd.read_csv(params.treatment_data, sep='\t', index_col=0)[['treatmentstatus', 'subtype_jpHMM', 'patientindex', 'sampledate_my']]
    treatment_df.index = treatment_df.index.map(str)
    treatment_df['sampledate_my'] = pd.to_datetime(treatment_df['sampledate_my'], infer_datetime_format=True)
    treatment_df = treatment_df.loc[~pd.isna(treatment_df['sampledate_my']), :]
    treatment_df['sampledate_my'] = treatment_df['sampledate_my'].apply(datetime2numeric)
    treatment_df = treatment_df.loc[treatment_df['subtype_jpHMM'] == params.subtype, ['treatmentstatus', 'patientindex', 'sampledate_my']]
    sample2patient = treatment_df.to_dict(orient='dict')['patientindex']

    id2tip = {_.name: _ for _ in tree}
    ids = sorted(list(treatment_df.index),
                 key=lambda _: (treatment_df.loc[_, 'patientindex'], treatment_df.loc[_, 'sampledate_my'],
                                _ not in id2tip))
    treatment_df = treatment_df.loc[ids, :]
    treatment_ids = set(treatment_df.index)

    arv_df = pd.read_csv(params.arv_data, sep='\t')
    arv_df['prefix'] = arv_df['mutation'].apply(lambda _: split_drm(_)[0])
    arv_df['letters'] = arv_df['mutation'].apply(lambda _: split_drm(_)[1])

    result_df = pd.DataFrame(columns=['class', 'reversion upper bounds', 'reversion lower bounds', 'reversion intervals',
                                      'reversion duration'])
    for drm in drms_by_pos:
        mp_df = pd.read_csv(params.mp_pattern.format(drm, drm), sep='\t', index_col=0)
        mp_df.index = mp_df.index.map(str)
        mp_ids = set(mp_df.index)

        prefix, letters = split_drm(drm)
        for letter in letters:
            resistant = '{}{}'.format(prefix, letter)
            print('----------', resistant, '----------')
            if 'polymorphic' == poly_df.loc[resistant, 'type']:
                print('skipping as polymorphic')
                continue

            # We need to check if a tree node id is in mp_ids because the upper part of the tree was cut as sensitive
            # and is not present in mp_df
            sensitive_node_ids = set(mp_df[(mp_df[resistant] <= 0.95) & (np.any(mp_df[:] > 0.95))].index) | (tree_ids - mp_ids)
            resistant_node_ids = set(mp_df[mp_df[resistant] > 0.95].index)

            id2date_last_sensitive = {}
            for t in tree:
                if t.name in resistant_node_ids:
                    n = t.up
                    while n.name not in sensitive_node_ids:
                        n = n.up
                    id2date_last_sensitive[sample2patient[t.name]] = getattr(n, DATE)
            rd_tree = set(id2date_last_sensitive.keys())

            resistant_sequence_ids = set(drm_df.loc[drm_df[resistant] == 'resistant', :].index) & treatment_ids
            sensitive_sequence_ids = set(drm_df.loc[drm_df[resistant] == 'sensitive', :].index) & treatment_ids

            id2date_first_resistant_naive = {}
            id2date_last_resistant_naive = {}
            id2date_first_sensitive = {}

            for patient in list(treatment_df.loc[list(resistant_sequence_ids), 'patientindex'].unique()):
                sample_ids = list(treatment_df.loc[treatment_df['patientindex'] == patient, :].index)
                first_res = None
                first_res_naive = None
                for _ in sample_ids:
                    sampling_date = treatment_df.loc[_, 'sampledate_my']
                    treatment_status = treatment_df.loc[_, 'treatmentstatus']
                    if _ in resistant_sequence_ids:
                        if first_res is None:
                            first_res = _
                        if 'naive' == treatment_status:
                            if first_res_naive is None:
                                first_res_naive = _
                                id2date_first_resistant_naive[patient] = sampling_date
                            else:
                                id2date_last_resistant_naive[patient] = sampling_date
                    else:
                        if first_res is None:
                            if _ in sensitive_sequence_ids:
                                if patient in rd_tree:
                                    print(id2date_last_sensitive[patient], '->',  sampling_date)
                                id2date_last_sensitive[patient] = sampling_date
                            continue
                        elif _ in sensitive_sequence_ids:
                            id2date_first_sensitive[patient] = sampling_date

            loss_durations = []
            censorship = []
            processed = set()
            for patient, date_last in id2date_last_resistant_naive.items():
                date_first = id2date_first_resistant_naive[patient]
                duration_kept = date_last - date_first
                if duration_kept:
                    processed.add(patient)
                    if patient in id2date_first_sensitive and patient in id2date_last_sensitive:
                        date_reverted = id2date_first_sensitive[patient]
                        date_acquired = id2date_last_sensitive[patient]
                        duration_lost = date_reverted - date_acquired
                        loss_durations.append([duration_kept, duration_lost])
                        # print(patient, date_acquired, date_first, date_last, date_reverted)
                        censorship.append(2)
                    else:
                        loss_durations.append(duration_kept)
                        censorship.append(1)
                        # print(patient, '--', date_first, date_last, '--')
            for patient in set(id2date_last_sensitive) - processed:
                if patient in id2date_first_sensitive:
                    date_reverted = id2date_first_sensitive[patient]
                    date_acquired = id2date_last_sensitive[patient]
                    duration_lost = date_reverted - date_acquired
                    loss_durations.append(duration_lost)
                    # print(patient, date_acquired, '-', '-', date_reverted)
                    censorship.append(-1)

            loss_durations = np.array(loss_durations, dtype=object)
            censorship = np.array(censorship)

            n_left = sum((censorship == -1).astype(int))
            n_right = sum((censorship == 1).astype(int))
            n_interval = sum((censorship == 2).astype(int))
            print('Found {} left-, {} right- and {} interval-censored data points'.format(n_left, n_right, n_interval))

            est, ci_min, ci_max = estimate_loss_time_and_ci(loss_durations, censorship)
            print(est, ci_min, ci_max)

            if est:
                result_df.loc[resistant.replace('_', ':'), :] = \
                    [get_arv_class(arv_df, prefix, letter),
                     n_left, n_right, n_interval,
                     '{:.1f} ({:.1f}--{:.1f})'.format(est, ci_min, ci_max) if est else '', ]

    result_df.join(pd.read_csv(params.Castro_data, sep='\t', index_col=0)).\
        to_csv(params.output_tab, sep='\t', index_label='DRM')