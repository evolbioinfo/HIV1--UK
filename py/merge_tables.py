import pandas as pd
from pastml.annotation import preannotate_forest
from pastml.tree import read_forest
from pastml.acr import _serialize_predicted_states, annotate_dates

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tabs', nargs='+', type=str)
    parser.add_argument('--tree', required=True, type=str)
    parser.add_argument('--output_tab', required=True, type=str)
    params = parser.parse_args()

    tree = read_forest(params.tree)[0]
    annotate_dates([tree])
    columns = []
    for tab in params.input_tabs:
        tab_df = pd.read_csv(tab, sep='\t', header=0, index_col=0)
        columns.extend(preannotate_forest(df=tab_df, forest=[tree])[0])
    _serialize_predicted_states(columns, params.output_tab, [tree])
    df = pd.read_csv(params.output_tab, sep='\t', index_col=0)
    df['consensus'] = ['+'.join(row[row != 'sensitive']) if sum((row != 'sensitive').astype(int)) > 0 else 'sensitive'
                       for (_, row) in df.iterrows()]
    df.to_csv(params.output_tab, sep='\t')

