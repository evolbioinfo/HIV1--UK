import pandas as pd
from pastml.tree import read_forest

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', required=True, type=str)
    parser.add_argument('--output', required=True, type=str)
    parser.add_argument('--tree', required=False, type=str, default=None)
    params = parser.parse_args()

    df = pd.read_table(params.input, sep='\t', index_col=0)
    df.index = df.index.map(str)
    df = df.loc[[_.name for _ in read_forest(params.tree)[0]], \
        [col for col in df.columns if col.startswith('RT') or col.startswith('PR') or col.startswith('IN')]]
    prevalence_df = pd.DataFrame((df == 'resistant').astype(int).sum().sort_values())
    prevalence_df.columns = ['Number of resistant cases']
    prevalence_df['Percentage of resistant cases'] = (100 * prevalence_df['Number of resistant cases'] / len(df))\
        .apply(lambda _: '{:.1f}'.format(_))
    prevalence_df.to_csv(params.output, index=True, sep='\t')


