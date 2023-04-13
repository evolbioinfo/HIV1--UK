import pandas as pd

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_data', required=True, type=str)
    parser.add_argument('--replacement_data', required=True, type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.input_data, sep='\t', index_col=0)
    df_repl = pd.read_csv(params.replacement_data, sep='\t', index_col=0)
    if len(df.columns) == len(df_repl.columns):
        df.columns = df_repl.columns
    else:
        for col in set(df_repl.columns) - set(df.columns):
            df[col] = 0
    common_index = list(set(df.index) & set(df_repl.index))
    df.loc[common_index, df.columns] = df_repl.loc[common_index, df.columns]
    df.to_csv(params.input_data, sep='\t', index_label='node')
