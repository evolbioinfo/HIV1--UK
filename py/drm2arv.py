import logging
import os.path
import re

import pandas as pd
import wikipedia

from sierrapy import SierraClient

QUERY = '''
drugResistance {
    gene { name },
    drugScores {
        drugClass { name },
        drug { displayAbbr, fullName },
        score,
        text,
        partialScores {
            mutations {
                text,
            },
        },
    }
}
'''


def get_date(drug):
    def find_dates_in_section(summary):
        for date in re.findall(r'\s[12][901][8901]\d[^\d]', summary):
            index = summary.find(date)
            prev_sentence = summary[max(0, index - 100): index]
            if 'approv' in prev_sentence or 'sold ' in prev_sentence:
                yield date[1:-1]

    page = wikipedia.page(wikipedia.search('{} HIV'.format(drug))[0])
    dates = set(find_dates_in_section(page.summary))
    if not dates:
        for _ in ('History', 'history'):
            if _ in page.sections:
                dates = set(find_dates_in_section(page.section(_)))
                if dates:
                    break
    if not dates and drug == 'rilpivirine':
        dates = {2011}
    return min(dates) if dates else None


if '__main__' == __name__:
    import argparse
    parser = argparse.ArgumentParser(description="Extracts SDRM drug resistance information.")
    parser.add_argument('--drm', nargs='+', default=['RT_A98G'], type=str, help="SDRMs of interest")
    parser.add_argument('--output', default='/home/azhukova/projects/hiv/data/metadata/arv_metadata.tab',
                        type=str, help="output file in tab-delimited format.")
    params = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    if os.path.exists(params.drm[0]):
        drms = set()
        for drm_file in params.drm:
            with open(drm_file, 'r') as f:
                drms |= set(f.read().strip().split(' '))
    else:
        drms = params.drm
    drms = [_.replace('PR_', 'PR:').replace('RT_', 'RT:').replace('IN_', 'IN:') for _ in sorted(drms)]

    data = []

    for dr in SierraClient().mutations_analysis(drms, QUERY)['drugResistance']:
        gene = dr['gene']['name']
        drug_scores = dr['drugScores']
        for ds in drug_scores:
            text = ds['text']
            if 'Resistance' not in text:
                continue
            drug_class = ds['drugClass']['name']
            drug_abbr = ds['drug']['displayAbbr'].replace('/r', '')
            drug_name = ds['drug']['fullName'].replace('/r', '')
            score = ds['score']
            partial_scores = ds['partialScores']
            for ps in partial_scores:
                for _ in ps['mutations']:
                    drm = _['text']
                    data.append(['{}:{}'.format(gene, drm), drug_class, drug_name, drug_abbr, score, text])
    df = pd.DataFrame(data, columns=['mutation', 'drug class', 'drug full name', 'drug abbreviation', 'score', 'note'])
    for drug in df['drug full name'].unique():
        df.loc[df['drug full name'] == drug, 'year'] = get_date(drug)
    df['mutation'] = df['mutation'].apply(lambda _: _.replace('PR:', 'PR_').replace('RT:', 'RT_'))

    df.sort_values(by=['mutation'], inplace=True)
    df.to_csv(params.output, sep='\t', index=False)


