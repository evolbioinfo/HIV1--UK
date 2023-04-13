import re

import pandas as pd
from ete3 import TreeNode
from pastml.tree import read_forest, annotate_dates, DATE, DATE_CI

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


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tree', required=True, type=str)
    parser.add_argument('--output_forest', required=True, type=str)
    parser.add_argument('--arv_tab', required=True, type=str)
    parser.add_argument('--poly', required=True, type=str)
    parser.add_argument('--drm', required=False, type=str, default=None)
    params = parser.parse_args()

    tree = read_forest(params.input_tree)[0]
    annotate_dates([tree])
    root_date = getattr(tree, DATE)
    print('Root year is {} ({}).'.format(root_date, getattr(tree, DATE_CI)))

    arv_df = pd.read_csv(params.arv_tab, index_col=None, sep='\t')

    prefix, letters = split_drm(params.drm)
    sub_drms = ['{}{}'.format(prefix, l) for l in letters]
    poly_status = set(pd.read_csv(params.poly, sep='\t', index_col=0).loc[sub_drms, 'type'].unique())
    if len(poly_status) > 1:
        raise ValueError('Polymorphic status is inconsistent for the mutations {}'.format(params.drm))
    if 'polymorphic' == poly_status.pop():
        print('{} is polymorphic'.format(params.drm))
        tree.write(outfile=params.output_forest, format=3, format_root_node=True, features=[DATE, DATE_CI])
    else:
        arv_year = FIRST_ARV_DATE
        arv_drm = next(_ for _ in arv_df['mutation'].unique() if _.startswith(prefix))
        arv_df = arv_df.loc[arv_df['mutation'] == arv_drm, :]
        if 'year' in arv_df.columns:
            dated_arv_df = arv_df.loc[~pd.isna(arv_df['year']), 'year']
            if len(dated_arv_df):
                arv_year = float(arv_df.loc[~pd.isna(arv_df['year']), 'year'].min())
        print('First {}-provoking ARV was accepted in {}.'.format(params.drm, arv_year))

        nwks = []
        todo = [(tree, root_date)]
        while todo:
            node, date = todo.pop()
            if date < arv_year:
                for child in node.children:
                    todo.append((child, date + child.dist))
            else:
                fake_root = TreeNode(dist=0, name='sensitive')
                fake_root.add_feature(DATE, arv_year)
                fake_root.add_feature(DATE_CI, [getattr(node.up, DATE), getattr(node, DATE)])
                fake_root.add_child(node, dist=date - arv_year)
                nwks.append(fake_root.write(format=3, format_root_node=True, features=[DATE, DATE_CI]))
        with open(params.output_forest, 'w+') as f:
            f.write('\n'.join(nwks))
