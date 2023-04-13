import os
import re

import pandas as pd
from pastml.annotation import preannotate_forest
from pastml.tree import read_forest, annotate_dates

def format_float_or_int(n):
    return (('{:.1f}' if n != int(n) else '{:.0f}').format(n)) if n else ''


def split_drm(drm):
    """
    Splits a given DRM into (prefix, to_letter)
    :param drm: input DRM
    :return: tuple (prefix, to_letter)
    """
    res = re.findall(r'((RT|PR|IN)[:_][A-Z][0-9]+)([A-Z]+)', drm)
    return (res[0][0], res[0][-1]) if res else (None, None)


def latexify(df):
    for col in ('class', '1st ARV'):
        df[col] = '\\scriptsize{' + df[col] + '}'
    df.index = '\\scriptsize{' + df.index + '}'
    for col in df.columns:
        df[col] = df[col].apply(lambda _: str(_).replace('%', '\\%').replace('(', '\\scriptsize{(').replace(')', ')}') \
            if not pd.isna(_) else '')


def format_year(year):
    return '\'{}'.format('{:.0f}'.format(year)[2:]) if year else ''


def get_drm_info(arv_df, prefix, drm, poly_df):
    drm_date, drug_class, drug = '', '', ''
    poly_status = poly_df.loc[drm, 'type']
    letter = drm.replace(prefix, '')
    sub_arv_df = arv_df[(arv_df['prefix'] == prefix) & arv_df['letters'].apply(lambda _: letter in _)]
    if 'year' in sub_arv_df.columns:
        dated_arv_df = sub_arv_df.loc[~pd.isna(sub_arv_df['year']), ['year', 'drug abbreviation', 'drug class']]
        if len(dated_arv_df):
            drm_date = float(dated_arv_df['year'].min())
            drug_class = dated_arv_df.loc[arv_df['year'] == drm_date, 'drug class'].min()
            drug = dated_arv_df.loc[arv_df['year'] == drm_date, 'drug abbreviation'].min()
    return drug_class, '{}{}'.format(drug, format_year(drm_date)) if 'polymorphic' != poly_status else 'polym.'


def count_tip_status(n):
    n_naive, n_treated, n_unknown = 0, 0, 0
    for t in n:
        status = getattr(t, 'treatmentstatus', None)
        if {'naive'} == status:
            n_naive += 1
        elif {'experienced'} == status:
            n_treated += 1
        else:
            n_unknown += 1
    return n_naive, n_treated, n_unknown


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--drms', type=str)
    parser.add_argument('--poly', type=str)
    parser.add_argument('--treatment_status', type=str)
    parser.add_argument('--arv_tab', type=str)
    parser.add_argument('--nwk', type=str)
    parser.add_argument('--output_tab', type=str)
    parser.add_argument('--mp_pattern', type=str)
    params = parser.parse_args()

    tree = read_forest(params.nwk, columns=params.drms)[0]
    annotate_dates([tree])
    treatment_df = pd.read_csv(params.treatment_status, sep='\t', index_col=0)[['treatmentstatus']]
    treatment_df.index = treatment_df.index.map(str)
    preannotate_forest([tree], treatment_df)
    poly_df = pd.read_csv(params.poly, sep='\t', index_col=0)
    arv_df = pd.read_csv(params.arv_tab, sep='\t')
    arv_df['prefix'] = arv_df['mutation'].apply(lambda _: split_drm(_)[0])
    arv_df['letters'] = arv_df['mutation'].apply(lambda _: split_drm(_)[1])

    result_df = pd.DataFrame(columns=['class', '1st ARV', 'n_res', 'resistant',
                                      'resistant-treated', 'resistant-naive',
                                      'TDR', 'TDR clusters', 'cluster sizes',
                                      'ADR', 'losses'])
    n_all = len(tree)

    with open(params.drms, 'r') as f:
        drms_by_pos = sorted(f.read().strip().split(' '))

    for drm in drms_by_pos:
        mp_df = pd.read_csv(params.mp_pattern.format(drm, drm), sep='\t', index_col=0)
        mp_df.index = mp_df.index.map(str)
        mp_ids = set(mp_df.index)

        prefix, letters = split_drm(drm)
        for letter in letters:
            resistant = '{}{}'.format(prefix, letter)
            is_polymorphic = poly_df.loc[resistant, 'type'] == 'polymorphic'

            resistant_node_ids = set(mp_df[mp_df[resistant] >= 0.5].index)

            resistant_tips = [_ for _ in tree if _.name in resistant_node_ids]
            n_res = len(resistant_tips)
            n_res_naive = sum(1 for _ in resistant_tips if getattr(_, 'treatmentstatus', None) == {'naive'})
            n_res_treated = sum(1 for _ in resistant_tips if getattr(_, 'treatmentstatus', None) == {'experienced'})

            n_adr, n_tdr, n_losses = 0, 0, 0
            n_singleton_clusters = 0
            tdr_roots = []
            for n in tree.traverse():
                # sensitive
                if n.name not in resistant_node_ids:
                    if not n.is_root() and n.up.name in resistant_node_ids:
                        n_losses += 1
                # resistant
                else:
                    if n.is_root() or n.up.name not in resistant_node_ids:
                        n_naive, n_treated, n_unknown = count_tip_status(n)
                        # probability that this cluster is naive-only
                        prob_all_naive = 0 if n_treated else 0.5 ** n_unknown
                        # if the mutation is polymorphic, we do not care whether the cluster is all-naive or not,
                        # as it could have as well appeared in a treatment-aive individual.
                        # So we will pretend that the source of resistance
                        # (the treated individual for non-polymorphic mutations) is identified.
                        if is_polymorphic:
                            prob_all_naive = 0
                        # hidden TDR
                        n_tdr += prob_all_naive
                        # ADR in an observed patient
                        n_adr += 1 - prob_all_naive
                        if n.is_leaf():
                            n_singleton_clusters += prob_all_naive
                        else:
                            tdr_roots.append(n)
                    if not n.is_leaf():
                        # In a (non-binary) tree,
                        # each internal node represents as many transmissions as it has children minus 1
                        n_tdr += len(n.children) - 1

            assert (n_adr + n_tdr - n_losses == n_res)

            cluster_sizes = []
            for r in tdr_roots:
                todo = [r]
                cs = 0
                while todo:
                    n = todo.pop()
                    if n.name not in resistant_node_ids:
                        continue
                    if n.is_leaf():
                        cs += 1
                    todo.extend(n.children)
                cluster_sizes.append(cs)

            n_clusters = len(tdr_roots) + n_singleton_clusters

            result_df.loc[resistant.replace('_', ':'), :] = \
                [*get_drm_info(arv_df, prefix=prefix, drm=resistant, poly_df=poly_df),
                 n_res,
                 '{r} ({p_r:.1f}%)'.format(r=n_res, p_r=100 * n_res / n_all),
                 '{r} ({p_r:.1f}%)'.format(r=n_res_treated, p_r=100 * n_res_treated / n_res) if n_res_treated else '',
                 '{r} ({p_r:.1f}%)'.format(r=n_res_naive, p_r=100 * n_res_naive / n_res) if n_res_naive else '',
                 '{tdr:.2f} ({p_tdr:.1f}%)'.format(tdr=n_tdr, p_tdr=100. * n_tdr / n_res) if n_tdr else '',
                 format_float_or_int(n_clusters),
                 '{cluster_min}-{cluster_max}'.format(cluster_min=1 if n_singleton_clusters else min(cluster_sizes),
                                                      cluster_max=max(cluster_sizes) if cluster_sizes else 1) \
                     if n_clusters else '',
                 '{adr:.2f} ({p_adr:.1f}%)'.format(adr=n_adr, p_adr=100. * n_adr / n_res) if n_adr else '',
                 '{loss} ({p_loss:.1f}%)'.format(loss=n_losses, p_loss=100. * n_losses / n_res) if n_losses else '',
                 ]
    result_df.sort_values(by=['n_res'], ascending=False, inplace=True)

    result_df.loc['RT:S68G', 'class'] = 'NRTI'
    latexify(result_df)
    result_df[['class', '1st ARV', 'resistant', 'resistant-treated', 'resistant-naive', 'TDR',
               'TDR clusters', 'cluster sizes', "ADR", 'losses']] \
        .to_csv(params.output_tab, sep='&', index_label='DRM')
    with open(params.output_tab, 'r') as f:
        lines = f.read().split('\n')
    with open(params.output_tab, 'w+') as f:
        f.write('\\\\\n'.join(lines))
