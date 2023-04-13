import gzip
import logging
import lzma
import os
import re
import pandas as pd
from Bio import SeqIO


def get_write_handle(fa, temp_suffix='', mode='wb'):
    if fa.endswith('.gz'):
        return gzip.open(fa + temp_suffix, mode)
    if fa.endswith('.xz'):
        return lzma.open(fa + temp_suffix, mode)
    return open(fa + temp_suffix, mode)


def get_read_handle(filepath, write=False):
    """Automatically get (un)compressed file handle"""
    ext = filepath.split(".")[-1]

    mode = "wt" if write else "rt"

    if ext == "gz":
        return gzip.open(filepath, mode)
    elif ext == "xz":
        return lzma.open(filepath, mode)

    return open(filepath, mode[0])


def clean_sequences(sequences, pos):
    for sequence in sequences:
        seq = str(sequence.seq)
        res = []
        end = 0
        for start in pos:
            res += seq[end: (start - 1) * 3]
            end = start * 3
        res += seq[end:]
        yield sequence.id, ''.join(res)


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_data', required=True, type=str)
    parser.add_argument('--input_fa', required=True, type=str)
    parser.add_argument('--output_fa', required=True, type=str)
    parser.add_argument('--PR_start_pos', required=True, type=int)
    parser.add_argument('--RT_start_pos', required=True, type=int)
    params = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    DRMS = list(pd.read_table(params.input_data, index_col=0).columns)
    pos = sorted({params.PR_start_pos + int(re.findall('\d+', _)[0]) for _ in DRMS if _.startswith('PR')}) \
          + sorted({params.RT_start_pos + int(re.findall('\d+', _)[0]) for _ in DRMS if _.startswith('RT')})

    with get_write_handle(params.output_fa, '.temp') as f:
        count = 0
        with get_read_handle(params.input_fa) as handle:
            for (seq_id, seq) in clean_sequences(SeqIO.parse(handle, "fasta"), pos):
                f.write('>{}\n{}\n'.format(seq_id, seq).encode())
                count += 1
    os.rename(params.output_fa + '.temp', params.output_fa)
    logging.info("Converted %d records to fasta" % count)



