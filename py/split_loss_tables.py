import re
import pandas as pd


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--B', type=str)
    parser.add_argument('--C', type=str)
    parser.add_argument('--tab', type=str)
    parser.add_argument('--sup_tab', type=str)
    params = parser.parse_args()

    df_B = pd.read_csv(params.B, sep='\t', index_col=[0, 1])
    df_C = pd.read_csv(params.C, sep='\t', index_col=[0, 1])
    df_C.drop(axis=1, labels=['DRM Castro et al.'], inplace=True)
    df = df_B.join(df_C, how='outer', lsuffix='B', rsuffix='C')
    for c in ('reversion upper bounds', 'reversion lower bounds', 'reversion intervals'):
        df[c + 'C'] = df[c + 'C'].apply(lambda _: '' if pd.isna(_) else int(_))
    df['reversion durationC'] = df['reversion durationC'].apply(lambda _: '' if pd.isna(_) else _)

    def key(v):
        enzyme, pos = re.findall(r'(PR|RT):[A-Z](\d+)', v[0])[0]
        return enzyme, int(pos)

    df['position'] = df.index.map(key)
    df.sort_values(by=['position'], inplace=True)

    df[['reversion durationB', 'reversion durationC', 'DRM Castro et al.']].to_csv(params.tab, sep='\t')
    df.drop(axis=1, labels=['reversion durationB', 'reversion durationC', 'DRM Castro et al.', 'position'], inplace=True)
    df.to_csv(params.sup_tab, sep='\t')

