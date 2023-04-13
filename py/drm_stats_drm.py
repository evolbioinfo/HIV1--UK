import os
import re

import pandas as pd
from pastml import numeric2datetime
from pastml.annotation import preannotate_forest
from pastml.tree import read_forest, annotate_dates, DATE, remove_certain_leaves

TABLE_HEADER = """\\begin{{table}}[!h] 
\\begin{{adjustwidth}}{{-\\extralength}}{{0cm}}
\\caption{{DRMs with prevalence $>0.5\\%$ found in position {position} in {subtype} data set, 
and the evolution of their presence over time.\\label{{tab:{position}}}}}
\\begin{{tabularx}}{{1.24\\textwidth}}{{crc|rrr|rrr|r|r}}
\\toprule
\\textbf{{date}} & \\textbf{{total}} & \\textbf{{DRM}} & \\multicolumn{{3}}{{c|}}{{\\textbf{{resistant cases}}}} & \\multicolumn{{3}}{{c|}}{{\\textbf{{TDR}}}} & \\multicolumn{{1}}{{c|}}{{\\textbf{{ADR}}}} & \\multicolumn{{1}}{{c}}{{\\textbf{{loss}}}}\\\\
& \\scriptsize{{samples}} & & \\multicolumn{{1}}{{c}}{{\\scriptsize{{(\\% of all)}}}} & \\multicolumn{{2}}{{|c|}}{{\\scriptsize{{treatment-}}}} & \\multicolumn{{1}}{{c}}{{\\scriptsize{{cases}}}} &  \\multicolumn{{2}}{{|c|}}{{\\scriptsize{{cluster}}}} & \\multicolumn{{1}}{{c|}}{{\\scriptsize{{cases}}}} & \\multicolumn{{1}}{{c}}{{\\scriptsize{{cases}}}}\\\\
& &  &  & \\multicolumn{{1}}{{|c}}{{\\scriptsize{{experienced}}}} & \\multicolumn{{1}}{{c}}{{\\scriptsize{{naive}}}}  & \\multicolumn{{1}}{{c}}{{\\scriptsize{{(\\% of}}}} &  \\multicolumn{{1}}{{|c}}{{\\scriptsize{{num.}}}} &  \\multicolumn{{1}}{{c|}}{{\\scriptsize{{sizes}}}} & \\multicolumn{{1}}{{c|}}{{\\scriptsize{{(\\% of}}}} & \\multicolumn{{1}}{{c}}{{\\scriptsize{{(\\% of}}}}\\\\
& & &  & \\multicolumn{{2}}{{|c|}}{{\\scriptsize{{(\\% of resistant)}}}} & \\multicolumn{{1}}{{c}}{{\\scriptsize{{resistant)}}}}  & \\multicolumn{{2}}{{|c|}}{{}} & \\multicolumn{{1}}{{c|}}{{\\scriptsize{{resistant)}}}} & \\multicolumn{{1}}{{c}}{{\\scriptsize{{resistant)}}}}\\\\
"""

TABLE_FOOTER = """\\bottomrule
\\end{tabularx}
\\end{adjustwidth}
\\end{table}
"""


def format_float_or_int(n):
    return (('{:.1f}' if n != int(n) else '{:.0f}').format(n)) if n else ''


def latexify(df, n_drms):
    col = 'year'
    df[col] = df[col].apply(lambda _: '\\midrule\\multirow{:.0f}'.format(n_drms) + '{*}{' + str(_) + '}' if _ and not pd.isna(_) else '')
    col = 'total'
    df[col] = df[col].apply(lambda _: '\\multirow{:.0f}'.format(n_drms) + '{*}{' + str(_) + '}' if _ and not pd.isna(_) else '')
    for col in df.columns:
        df[col] = df[col].apply(lambda _: str(_).replace('%', '\\%').replace('(', '\\scriptsize{(').replace(')', ')}') \
            if not pd.isna(_) else '')


def count_tip_status(n):
    n_naive, n_treated, n_unknown = 0, 0, 0
    for t in n:
        status = getattr(t, 'treatmentstatus', None)
        if {'naive'} == status:
            n_naive += 1
        elif {'experienced'} == status:
            n_treated += 1
        else:
            n_unknown += 1
    return n_naive, n_treated, n_unknown


def split_drm(drm):
    """
    Splits a given DRM into (prefix, to_letter)
    :param drm: input DRM
    :return: tuple (prefix, to_letter)
    """
    res = re.findall(r'((RT|PR|IN)[:_][A-Z][0-9]+)([A-Z]+)', drm)
    return (res[0][0], res[0][-1]) if res else (None, None)


def format_year(year):
    return '\'{}'.format('{:.0f}'.format(year)[2:]) if year else ''


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--treatment_status', type=str)
    parser.add_argument('--nwk', type=str)
    parser.add_argument('--output_tab', type=str)
    parser.add_argument('--mp', type=str)
    parser.add_argument('--drm', type=str)
    parser.add_argument('--years', type=int, nargs='+')
    parser.add_argument('--subtype', type=str)
    params = parser.parse_args()

    tree = read_forest(params.nwk)[0]
    annotate_dates([tree])
    treatment_df = pd.read_csv(params.treatment_status, sep='\t', index_col=0)[['treatmentstatus']]
    treatment_df.index = treatment_df.index.map(str)
    preannotate_forest([tree], treatment_df)

    result_df = pd.DataFrame(columns=['year', 'total', 'DRM', 'resistant',
                                      'resistant-treated', 'resistant-naive',
                                      'TDR', 'TDR clusters', 'cluster sizes',
                                      'ADR', 'losses'])

    mp_df = pd.read_csv(params.mp, sep='\t', index_col=0)
    mp_df.index = mp_df.index.map(str)
    mp_ids = set(mp_df.index)

    position, letters = split_drm(params.drm)

    letter2resistant_node_ids = {}
    for letter in letters:
        resistant = '{}{}'.format(position, letter)
        letter2resistant_node_ids[letter] = set(mp_df[mp_df[resistant] > 0.5].index)

    for year in reversed(sorted(params.years)):
        tree = remove_certain_leaves(tree, lambda _: getattr(_, DATE) > year)
        if not tree:
            continue
        print(year, len(tree))
        date = numeric2datetime(max(getattr(_, DATE) for _ in tree)).strftime("%d-%m-%y")

        tips = list(tree)
        n_tips = len(tips)
        for letter in letters:
            resistant_node_ids = letter2resistant_node_ids[letter]

            resistant_tips = [_ for _ in tree if _.name in resistant_node_ids]
            n_res = len(resistant_tips)
            n_res_naive = sum(1 for _ in resistant_tips if getattr(_, 'treatmentstatus', None) == {'naive'})
            n_res_treated = sum(1 for _ in resistant_tips if getattr(_, 'treatmentstatus', None) == {'experienced'})

            n_adr, n_tdr, n_losses = 0, 0, 0
            n_singleton_clusters = 0
            tdr_roots = []
            for n in tree.traverse():
                # sensitive
                if n.name not in resistant_node_ids:
                    if not n.is_root() and n.up.name in resistant_node_ids:
                        n_losses += 1
                # resistant
                else:
                    if n.is_root() or n.up.name not in resistant_node_ids:
                        n_naive, n_treated, n_unknown = count_tip_status(n)
                        # probability that this cluster is naive-only
                        prob_all_naive = 0 if n_treated else 0.5 ** n_unknown
                        # hidden TDR
                        n_tdr += prob_all_naive
                        # ADR in an observed patient
                        n_adr += 1 - prob_all_naive
                        if n.is_leaf():
                            n_singleton_clusters += prob_all_naive
                        else:
                            tdr_roots.append(n)
                    if not n.is_leaf():
                        # In a (non-binary) tree,
                        # each internal node represents as many transmissions as it has children minus 1
                        n_tdr += len(n.children) - 1

            assert (n_adr + n_tdr - n_losses == n_res)

            cluster_sizes = []
            for r in tdr_roots:
                todo = [r]
                cs = 0
                while todo:
                    n = todo.pop()
                    if n.name not in resistant_node_ids:
                        continue
                    if n.is_leaf():
                        cs += 1
                    todo.extend(n.children)
                cluster_sizes.append(cs)

            n_clusters = len(tdr_roots) + n_singleton_clusters

            result_df.loc[date+letter, :] = \
                [date if letter == letters[0] else '',
                 n_tips if letter == letters[0] else '',
                 letter,
                 '{r} ({p_r:.1f}%)'.format(r=n_res, p_r=100 * n_res / n_tips) if n_res else '',
                 '{r} ({p_r:.1f}%)'.format(r=n_res_treated, p_r=100 * n_res_treated / n_res) if n_res_treated else '',
                 '{r} ({p_r:.1f}%)'.format(r=n_res_naive, p_r=100 * n_res_naive / n_res) if n_res_naive else '',
                 '{tdr:.2f} ({p_tdr:.1f}%)'.format(tdr=n_tdr, p_tdr=100. * n_tdr / n_res) if n_tdr else '',
                 format_float_or_int(n_clusters),
                 '{cluster_min}-{cluster_max}'.format(cluster_min=1 if n_singleton_clusters else min(cluster_sizes),
                                                      cluster_max=max(cluster_sizes) if cluster_sizes else 1) \
                     if n_clusters else '',
                 '{adr:.2f} ({p_adr:.1f}%)'.format(adr=n_adr, p_adr=100. * n_adr / n_res) if n_adr else '',
                 '{loss} ({p_loss:.1f}%)'.format(loss=n_losses, p_loss=100. * n_losses / n_res) if n_losses else '']
    latexify(result_df, len(letters))
    result_df.to_csv(params.output_tab, sep='&', index=False, header=False)

    with open(params.output_tab, 'r') as f:
        lines = f.read().split('\n')
    with open(params.output_tab, 'w+') as f:
        f.write(TABLE_HEADER.format(position=position.replace('_', ':'), subtype=params.subtype))
        f.write('\\\\\n'.join(lines))
        f.write(TABLE_FOOTER)



