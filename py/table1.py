import re

import pandas as pd
from pastml.tree import read_tree

POLYMORPHIC = 'polymorphic'

NONPOLYMORPHIC = 'nonpolymorphic'

SUBTYPE = 'subtype_jpHMM'

PATIENTINDEX = 'patientindex'

TREATMENTSTATUS = 'treatmentstatus'


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--metadata', default='/home/azhukova/projects/hiv/data/metadata/metadata.uk.tab', type=str)
    parser.add_argument('--metadata_drm', default=['/home/azhukova/projects/hiv/data/metadata/B/metadata.drm.tab',
                                               '/home/azhukova/projects/hiv/data/metadata/C/metadata.drm.tab'],
                        nargs='+', type=str)
    parser.add_argument('--poly', default='/home/azhukova/projects/hiv/data/metadata/drm_types.tab', type=str)
    parser.add_argument('--output', default='/home/azhukova/projects/hiv/data/metadata/Table1.tab', type=str)
    parser.add_argument('--tree', default=['/home/azhukova/projects/hiv/data/trees/B/rooted_raxml.nwk',
                                           '/home/azhukova/projects/hiv/data/trees/C/rooted_raxml.nwk'],
                        nargs='+', type=str)
    parser.add_argument('--timetree', default=['/home/azhukova/projects/hiv/data/timetrees/B/raxml.lsd2.nwk',
                                               '/home/azhukova/projects/hiv/data/timetrees/C/raxml.lsd2.nwk'],
                        nargs='+', type=str)
    parser.add_argument('--subtype', nargs='+', default=['B', 'C'], type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.metadata, sep='\t', index_col=0)[[TREATMENTSTATUS, PATIENTINDEX, SUBTYPE]]
    df.index = df.index.map(str)
    results = pd.DataFrame(columns=params.subtype)
    
    poly = pd.read_csv(params.poly, sep='\t', index_col=0)
    non_poly_drms = set(poly[poly['type'] == NONPOLYMORPHIC].index)
    poly_drms = set(poly[poly['type'] == POLYMORPHIC].index)
    
    for (subtype, nwk, dated_nwk, md) in zip(params.subtype, params.tree, params.timetree, params.metadata_drm):
        st_df = df[df[SUBTYPE] == subtype]
        n_total = len(st_df)
        results.loc['total', subtype] = n_total
        tree = read_tree(nwk)
        # tree length divided by the number of branches
        phyl_diversity = sum(_.dist for _ in tree.traverse()) / (2 * len(tree) - 3)
        n_firstseq = len(tree)
        results.loc['filtered by patient (first only, % of total)', subtype] = \
            '{} ({:.0f}%)'.format(n_firstseq, 100 * n_firstseq / n_total)
        tree = read_tree(dated_nwk)
        n_wo_outliers = len(tree)
        results.loc['without temporal outliers (% of filtered)', subtype] = \
            '{} ({:.0f}%)'.format(n_wo_outliers, 100 * n_wo_outliers / n_firstseq)
    
        tdf = st_df.loc[[_.name for _ in tree], [TREATMENTSTATUS]]
    
        drm_df = pd.read_csv(md, sep='\t', index_col=0)
        drm_df.index = drm_df.index.map(str)
        drm_df = drm_df.loc[[_.name for _ in tree], \
            [col for col in drm_df.columns if 'PR' in col or 'RT' in col or 'IN' in col]]
        nonpoly_df = drm_df.loc[:, [c for c in drm_df.columns if c in non_poly_drms]]
        poly_df = drm_df.loc[:, [c for c in drm_df.columns if c in poly_drms]]
    
        resistant_count_df = pd.DataFrame((drm_df == 'resistant').astype(int).sum(axis=1))
        resistant_count_df.columns = ['Number of resistant mutations']
        resistant_count_df = resistant_count_df[resistant_count_df['Number of resistant mutations'] > 0]
        n_resistant = len(resistant_count_df)
        results.loc['with DRM(s) (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n_resistant, 100 * n_resistant / n_wo_outliers)
        n1drm = len(resistant_count_df[resistant_count_df['Number of resistant mutations'] == 1])
        results.loc['w. 1 DRM (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n1drm, 100 * n1drm / n_wo_outliers)
        n2plus_drm = len(resistant_count_df[resistant_count_df['Number of resistant mutations'] > 1])
        results.loc['w. >=2 DRMs (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n2plus_drm, 100 * n2plus_drm / n_wo_outliers)
        res_ids = set(resistant_count_df.index)
    
        nonpoly_count_df = pd.DataFrame((nonpoly_df == 'resistant').astype(int).sum(axis=1))
        nonpoly_count_df.columns = ['Number of resistant mutations']
        nonpoly_count_df = nonpoly_count_df[nonpoly_count_df['Number of resistant mutations'] > 0]
        n_nonpoly = len(nonpoly_count_df)
        results.loc['with np DRM(s) (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n_nonpoly, 100 * n_nonpoly / n_wo_outliers)
        n1drm = len(nonpoly_count_df[nonpoly_count_df['Number of resistant mutations'] == 1])
        results.loc['w. 1 np DRM (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n1drm, 100 * n1drm / n_wo_outliers)
        n2plus_drm = len(nonpoly_count_df[nonpoly_count_df['Number of resistant mutations'] > 1])
        results.loc['w. >=2 np DRMs (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n2plus_drm, 100 * n2plus_drm / n_wo_outliers)
        nonpoly_ids = set(nonpoly_count_df.index)
    
        poly_count_df = pd.DataFrame((poly_df == 'resistant').astype(int).sum(axis=1))
        poly_count_df.columns = ['Number of resistant mutations']
        poly_count_df = poly_count_df[poly_count_df['Number of resistant mutations'] > 0]
        n_poly = len(poly_count_df)
        results.loc['with p DRM(s) (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n_poly, 100 * n_poly / n_wo_outliers)
        n1drm = len(poly_count_df[poly_count_df['Number of resistant mutations'] == 1])
        results.loc['w. 1 p DRM (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n1drm, 100 * n1drm / n_wo_outliers)
        n2plus_drm = len(poly_count_df[poly_count_df['Number of resistant mutations'] > 1])
        results.loc['w. >=2 p DRMs (% of w/o outliers)', subtype] = \
            '{} ({:.0f}%)'.format(n2plus_drm, 100 * n2plus_drm / n_wo_outliers)
        poly_ids = set(poly_count_df.index)
    
        for status in ('naive', 'experienced', 'unknown'):
            status_tdf = \
                tdf[(tdf[TREATMENTSTATUS] == status) if status != 'unknown' else pd.isna(tdf[TREATMENTSTATUS])]
            n_status = len(status_tdf)
            results.loc['treatment-{}  (% of w/o outliers)'.format(status), subtype] = \
                '{} ({:.0f}%)'.format(n_status, 100 * n_status / n_wo_outliers)
            if 'unknown' != status:
                n_status_res = len(set(status_tdf.index) & res_ids)
                results.loc['treatment-{} with DRM(s)  (% of tr.-{})'.format(status, status), subtype] = \
                    '{} ({:.0f}%)'.format(n_status_res, 100 * n_status_res / n_status)
                n_status_nonpoly = len(set(status_tdf.index) & nonpoly_ids)
                results.loc['treatment-{} with np DRM(s)  (% of tr.-{})'.format(status, status), subtype] = \
                    '{} ({:.0f}%)'.format(n_status_nonpoly, 100 * n_status_nonpoly / n_status)
                n_status_poly = len(set(status_tdf.index) & poly_ids)
                results.loc['treatment-{} with p DRM(s)  (% of tr.-{})'.format(status, status), subtype] = \
                    '{} ({:.0f}%)'.format(n_status_poly, 100 * n_status_poly / n_status)
                
        with open(dated_nwk.replace('.nwk', '.log'), 'r') as f:
            r, r_min, r_max, d, d_min, d_max = None, None, None, None, None, None
            for line in f:
                rate = re.findall(r'rate (\d+\.\d+) \[(\d+\.\d+); (\d+\.\d+)\], '
                                  r'tMRCA (\d\d\d\d)\.\d+ \[\d\d(\d\d)\.\d+; \d\d(\d\d)\.\d+\]', line)
                if rate:
                    r, r_min, r_max, d, d_min, d_max = rate[0]
        results.loc['Root date', subtype] = "{} ('{}-'{})".format(int(d), d_min, d_max)
        results.loc['Mutation rate x 10-3', subtype] = "{:.1f} ({:.1f}-{:.1f})".format(float(r) * 1000,
                                                                                         float(r_min) * 1000,
                                                                                         float(r_max) * 1000)
        results.loc['Phylogenetic diversity', subtype] = '{:.3f}'.format(phyl_diversity)
    
    results.to_csv(params.output,index=True, sep='\t')