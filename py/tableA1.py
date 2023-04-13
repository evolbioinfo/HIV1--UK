import numpy as np
import pandas as pd
from pastml.tree import read_tree

NONPOLYMORPHIC = 'nonpolymorphic'

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--metadata', default=['/home/azhukova/projects/hiv/data/metadata/B/metadata.drm.tab',
                                               '/home/azhukova/projects/hiv/data/metadata/C/metadata.drm.tab'],
                        nargs='+', type=str)
    parser.add_argument('--poly', default='/home/azhukova/projects/hiv/data/metadata/drm_types.tab', type=str)
    parser.add_argument('--output', default='/home/azhukova/projects/hiv/data/metadata/Table2.tab', type=str)
    parser.add_argument('--tree', default=['/home/azhukova/projects/hiv/data/timetrees/B/raxml.lsd2.nwk',
                                           '/home/azhukova/projects/hiv/data/timetrees/C/raxml.lsd2.nwk'],
                        nargs='+', type=str)
    parser.add_argument('--subtype', nargs='+', default=['B', 'C'], type=str)
    params = parser.parse_args()

    poly = pd.read_csv(params.poly, sep='\t', index_col=0)
    non_poly_drms = set(poly[poly['type'] == NONPOLYMORPHIC].index)
    result = None
    for (subtype, nwk, md) in zip(params.subtype, params.tree, params.metadata):
        df = pd.read_csv(md, sep='\t', index_col=0)
        df.index = df.index.map(str)
        tree = read_tree(nwk)
        df = df.loc[[_.name for _ in tree], [col for col in df.columns if 'PR' in col or 'RT' in col or 'IN' in col]]

        resistant_count_df = pd.DataFrame((df == 'resistant').astype(int).sum(axis=1))
        resistant_count_df.columns = ['Number of resistant mutations']
        resistant_count_df['Number of cases {}'.format(subtype)] = 1

        drm_count = resistant_count_df.groupby(by='Number of resistant mutations').count()
        if result is None:
            result = drm_count
        else:
            result = result.join(drm_count, how="outer")

        df = df.loc[:, [c for c in df.columns if c in non_poly_drms]]

        resistant_count_df = pd.DataFrame((df == 'resistant').astype(int).sum(axis=1))
        resistant_count_df.columns = ['Number of resistant mutations']
        resistant_count_df['Number of cases {} (nonpolymorphic only)'.format(subtype)] = 1

        result = result.join(resistant_count_df.groupby(by='Number of resistant mutations').count(), how="outer")

    counter_B, counter_C, counter_np_B, counter_np_C = [], [], [], []

    def get_num(val):
        if pd.isna(val):
            return 0
        return int(val)

    for num, row in result.iterrows():
        if num > 0:
            counter_B.extend([num] * get_num(row['Number of cases B']))
            counter_C.extend([num] * get_num(row['Number of cases C']))
            counter_np_B.extend([num] * get_num(row['Number of cases B (nonpolymorphic only)']))
            counter_np_C.extend([num] * get_num(row['Number of cases C (nonpolymorphic only)']))
    print(np.median(counter_B), np.median(counter_C), np.median(counter_np_B), np.median(counter_np_C))


    for col in result.columns:
        result[col].fillna(0, inplace=True)
        total = sum(result[col])
        result[col] = result[col].astype(int).apply(lambda _: '{} ({:.2f}%)'.format(_, 100 * _ / total))
    result.to_csv(params.output, index=True, sep='\t')

