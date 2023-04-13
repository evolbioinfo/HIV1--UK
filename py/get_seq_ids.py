import logging
import pandas as pd


PATIENT_COL = 'patientindex'
SAMPLEDATE_COL = 'sampledate_my'

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_data', required=True, type=str)
    parser.add_argument('--output_data', required=False, type=str)
    parser.add_argument('--outgroup', required=False, type=str)
    parser.add_argument('--subtype', required=True, type=str)
    parser.add_argument('--first_sample', action='store_true')
    params = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    df = pd.read_csv(params.input_data, sep='\t', index_col=0, header=0)
    df.index = df.index.map(str)

    if params.first_sample:
        df[SAMPLEDATE_COL] = pd.to_datetime(df[SAMPLEDATE_COL])
        df.sort_values(by=SAMPLEDATE_COL, ascending=True, inplace=True)
        old_len = len(df)
        n = len(df)
        df.drop_duplicates(subset=[PATIENT_COL], inplace=True)
        print('Kept {} first sequences (out of {})'.format(len(df), n))

    if params.output_data:
        in_df = df[df['compatible_subtypes'].apply(lambda _: (params.subtype in _.split('/')) if not pd.isna(_) else False)]
        logging.info('Extracted {} ids of subtype {}'.format(len(in_df), params.subtype))
        with open(params.output_data, 'w+') as f:
            f.write('\n'.join(list(in_df.index.map(str))))
    if params.outgroup:
        out_df = df[df['subtype_jpHMM'].apply(lambda _: (params.subtype != _ and len(_) == 1) if not pd.isna(_) else False)].sample(5)
        with open(params.outgroup, 'w+') as f:
            f.write('\n'.join(list(out_df.index.map(str))))
