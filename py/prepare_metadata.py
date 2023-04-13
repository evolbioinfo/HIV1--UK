import pandas as pd

TREATMENT_COL = 'treatmentstatus'
DATE_COL = 'sampledate_my'
PATIENT_COL = 'patientindex'

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_data', required=True, type=str)
    parser.add_argument('--output_data', required=True, type=str)
    parser.add_argument('--subtype_data', nargs='*', type=str, default=[])
    params = parser.parse_args()

    df = pd.read_csv(params.input_data, index_col=1, header=0)
    df.index = df.index.map(str)

    # rename treatment status
    cat2name = {1: 'naive', 2: 'experienced', 3: None}
    df[TREATMENT_COL] = df[TREATMENT_COL].map(lambda _: cat2name[_])
    df[[PATIENT_COL, DATE_COL, TREATMENT_COL]].to_csv(params.output_data, sep='\t', header=True, index=True, index_label='id')
