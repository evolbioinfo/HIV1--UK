import re
import shutil

import pandas as pd
from pastml.tree import read_forest, annotate_dates, DATE

# Date of AZT acceptance
FIRST_ARV_DATE = 1987


def split_drm(drm):
    """
    Splits a given DRM into (prefix, to_letter)
    :param drm: input DRM
    :return: tuple (prefix, to_letter)
    """
    res = re.findall(r'((RT|PR|IN)[:_][A-Z][0-9]+)([A-Z]+)', drm)
    return (res[0][0], res[0][-1]) if res else (None, None)


def get_drm_date(arv_tab, drm):
    arv_df = pd.read_csv(arv_tab, sep='\t')
    prefix, letters = split_drm(drm)
    arv_drm = next(_ for _ in arv_df['mutation'].unique() if _.startswith(prefix))
    arv_df = arv_df.loc[arv_df['mutation'] == arv_drm, :]
    drm_date = FIRST_ARV_DATE
    if 'year' in arv_df.columns:
        dated_arv_df = arv_df.loc[~pd.isna(arv_df['year']), 'year']
        if len(dated_arv_df):
            drm_date = float(arv_df.loc[~pd.isna(arv_df['year']), 'year'].min())
    return drm_date


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tab', required=True, type=str)
    parser.add_argument('--input_tree', required=True, type=str)
    parser.add_argument('--output_tab', required=True, type=str)
    parser.add_argument('--arv_tab', required=True, type=str)
    parser.add_argument('--poly', required=True, type=str)
    parser.add_argument('--drm', required=True, type=str)
    params = parser.parse_args()

    prefix, letters = split_drm(params.drm)
    sub_drms = ['{}{}'.format(prefix, l) for l in letters]
    poly_status = set(pd.read_csv(params.poly, sep='\t', index_col=0).loc[sub_drms, 'type'].unique())
    if len(poly_status) > 1:
        raise ValueError('Polymorphic status is inconsistent for the mutations {}'.format(params.drm))
    if 'polymorphic' == poly_status.pop():
        print('{} is polymorphic'.format(params.drm))
        shutil.copyfile(params.input_tab, params.output_tab)
    else:
        tree = read_forest(params.input_tree)[0]
        annotate_dates([tree])

        acr_df = pd.read_csv(params.input_tab, index_col=0, sep='\t')
        acr_df.index = acr_df.index.map(str)
        acr_df = acr_df.loc[acr_df.index != 'sensitive', :]

        drm_date = get_drm_date(params.arv_tab, params.drm)

        for n in tree.traverse():
            date = getattr(n, DATE)
            if date <= drm_date:
                acr_df.loc[n.name, params.drm] = 'sensitive'

        acr_df.to_csv(params.output_tab, sep='\t', index_label='node')

