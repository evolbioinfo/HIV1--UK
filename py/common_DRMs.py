import pandas as pd


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', required=True, type=str)
    parser.add_argument('--output', required=True, type=str)
    parser.add_argument('--threshold', default=0.5, type=float)
    params = parser.parse_args()

    df = pd.read_table(params.input, index_col=0)
    with open(params.output, 'w+') as f:
        f.write(' '.join(df[df['Percentage of resistant cases'] > params.threshold].index))

