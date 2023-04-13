import os

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--B', type=str, nargs='+')
    parser.add_argument('--C', type=str, nargs='+')
    parser.add_argument('--tex', type=str)
    parser.add_argument('--pattern', type=str)
    params = parser.parse_args()

    with open(params.tex, 'w+') as f:
        f.write('\\section{Resistance statistics over time by DRM}\n')
        for i, subtype in enumerate(['B', 'C']):
            if i:
                f.write('\\pagebreak\n')
            f.write('\\subsection{{Subtype {subtype}}}\n'.format(subtype=subtype))
            for drm in (params.B if 'B' == subtype else params.C):
                table_file = os.path.join(params.pattern.format(subtype, drm))
                with open(table_file, 'r') as tf:
                    f.write(tf.read())
                f.write('\n\n')
