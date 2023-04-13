import re
from collections import defaultdict

import pandas as pd


def split_drm(drm):
    """
    Splits a given DRM into (prefix, to_letter)
    :param drm: input DRM
    :return: tuple (prefix, to_letter)
    """
    res = re.findall(r'((RT|PR|IN)[:_][A-Z][0-9]+)([A-Z])', drm)
    return (res[0][0], res[0][-1]) if res else (None, None)


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_data', default='/home/azhukova/projects/hiv/data/metadata/C/metadata.drm.tab', type=str)
    parser.add_argument('--output_data', default='/home/azhukova/projects/hiv/data/metadata/C/metadata.drm.by_position.tab', type=str)
    parser.add_argument('--input_common_drms', default='/home/azhukova/projects/hiv/data/metadata/C/common_drms.txt', type=str)
    parser.add_argument('--output_common_drms', default='/home/azhukova/projects/hiv/data/metadata/C/common_drms.by_position.txt', type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.input_data, sep='\t', index_col=0)
    with open(params.input_common_drms, 'r') as f:
        columns = f.read().strip().split(' ')

    pos2letters = defaultdict(set)
    for col in columns:
        (prefix, to_letter) = split_drm(col)
        if prefix:
            pos2letters[prefix].add(to_letter)

    with open(params.output_common_drms, 'w+') as f:
        f.write(' '.join(sorted('{}{}'.format(k, ''.join(sorted(vs))) for (k, vs) in pos2letters.items())))

    result_df = pd.DataFrame(index=list(df.index))
    for col in sorted(columns):
        (prefix, to_letter) = split_drm(col)
        if prefix:
            letters = pos2letters[prefix]
            new_col = '{}{}'.format(prefix, ''.join(sorted(letters)))
            if new_col not in result_df.columns:
                result_df[new_col] = 'sensitive'
            result_df.loc[df[col] == 'resistant', new_col] = col
            result_df.loc[pd.isna(df[col]), new_col] = None

    result_df.index = result_df.index.map(str)
    result_df.loc['sensitive', :] = 'sensitive'

    result_df.to_csv(params.output_data, sep='\t', index_label='inputSequence')
