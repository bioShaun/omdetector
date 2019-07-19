import re
import fire
import numpy as np
import pandas as pd

BLAST_HEADER = [
    'qseqid',
    'sseqid',
    'pident',
    'qlen',
    'length',
    'mismatch',
    'gapopen',
    'evalue',
    'bitscore',
    'stitle',
]

SWISS_PATTERN = re.compile(r"OS=(.*?) \w+=")
UNIREF_PATTERN = re.compile(r"Tax=(.*?) \w+=")


def sp_name(stitle):
    swiss_title = SWISS_PATTERN.search(stitle)
    uniref_title = UNIREF_PATTERN.search(stitle)
    if swiss_title:
        return swiss_title.groups()[0]
    elif uniref_title:
        return uniref_title.groups()[0]
    else:
        return np.nan


def sp_portion(blastout, outfile):
    blast_df = pd.read_csv(blastout, header=None, names=BLAST_HEADER, sep='\t')
    target_des_df = blast_df.loc[:, ['qseqid', 'stitle']].drop_duplicates()
    target_des_df.loc[:, 'species'] = target_des_df.stitle.map(sp_name)
    target_des_df.dropna(inplace=True)
    target_des_df.loc[:, 'tr_id'] = target_des_df.qseqid.map(lambda x: ''.join(
        x.split('.')[:-1]))
    tr_sp_count = target_des_df.groupby(['tr_id'])['species'].value_counts()
    orf_count = tr_sp_count.sum(level=0)
    tr_sp_rate = tr_sp_count / orf_count
    tr_sp_rate.name = 'number'
    tr_sp_rate = tr_sp_rate.reset_index()
    sp_count = tr_sp_rate.groupby(['species'])['number'].sum()
    sp_count.sort_values(ascending=False, inplace=True)
    sp_portion = sp_count / sp_count.sum()
    sp_stats = pd.concat([sp_count, sp_portion], axis=1)
    sp_stats.columns = ['hit_count', 'hit_portion']
    sp_stats.to_csv(outfile, sep='\t')


if __name__ == "__main__":
    fire.Fire(sp_portion)
