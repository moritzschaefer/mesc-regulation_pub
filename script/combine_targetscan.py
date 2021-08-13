from functools import partial
from multiprocessing import Pool

import gffutils
import pandas as pd
from tqdm import tqdm


def cut_version(s: str):
    if type(s) is str:
        if '.' in s:
            return s[:s.find('.')]
        else:
            return s
    else:
        try:
            name = s.name
        except AttributeError:
            name = None
        return pd.Series([cut_version(v) for v in s], name=name)
conserved_df = pd.read_csv(snakemake.input['conserved'][0], sep='\t')
nonconserved_df = pd.read_csv(snakemake.input['nonconserved'][0], sep='\t')

# output_columns = [
#         'Gene Symbol', 'miRNA', 'context++ score',
#         'context++ score percentile', 'context++ score',
#         'context++ score percentile', 'geneId'
#     ]

conserved_df['conserved'] = 1
nonconserved_df['conserved'] = 0
df = pd.concat([conserved_df, nonconserved_df], ignore_index=True)

df['mre_type'] = df['Site Type'].apply({1: '7merA1', 2: '7merm8', 3: '8mer'}.get)
df['Geneid'] = cut_version(df['Gene ID'])
df.drop(df.index[df['Geneid'].isna()], inplace=True)
df.drop(df.index[df['Gene Tax ID'] != 10090], inplace=True)
# approximately 10% of nonconserved don't have a context++ score. delete them
df.drop(df.index[df['context++ score'].isna()], inplace=True)
df.drop(columns=['Gene Tax ID', 'Gene ID', 'Site Type'], inplace=True)

# find genome-position of MRE
db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)
df['mre_start'] = -1
df['mre_end'] = -1

grouped = df.groupby('Transcript ID')
def iter_mres(args):
    index, utr_start, utr_end, strand, threeputr_features = args
    # for index, row in group.iterrows():
    pos = 0
    for tpu_len, tpu_start, tpu_end in threeputr_features:
        if pos + tpu_len > utr_start:
            within_feature = utr_start - pos
            within_feature_end = utr_end - pos
            if strand == '+':
                return (index, tpu_start + within_feature, tpu_start + within_feature_end)
            else:
                return (index, tpu_end - within_feature_end, tpu_end - within_feature)
        pos += tpu_len
    # if the position of the MRE is *after* the annotated 3pUTR, then just imaginary extend the last exon.
    # HOWEVER we return the value multiplied by -1, to be able to count it later
    try:
        if strand == '+':
            return (index, tpu_end + (utr_start - pos), tpu_end + (utr_end - pos))
        else:
            return (index, tpu_start - (utr_end - pos), tpu_start - (utr_start - pos))
    except NameError:
        return (index, -1, -1)

def find_threeputr_features(transcript):
    '''
    should anyways only return one at most. Maybe I should assert this..
    '''
    threeputr_features = [(len(tpu), tpu.start, tpu.end) for tpu in
                          db.children(transcript, featuretype='three_prime_utr', order_by='start' if transcript.strand == '+' else 'end')]
    if len(threeputr_features) == 0:  # try to find UTR features
        threeputr_features = [(len(tpu), tpu.start, tpu.end) for tpu in
         db.children(transcript, featuretype='UTR', order_by='start' if transcript.strand == '+' else 'end')
         if (transcript.strand == '+' and tpu.end == transcript.end)  # only select 3'UTR.
         or (transcript.strand == '-' and tpu.start == transcript.start)
        ]
    return threeputr_features

mre_starts = {}
mre_ends = {}
with Pool(processes=snakemake.threads) as pool:
    for tid, group in tqdm(grouped, total=len(grouped)):
        try:
            transcript = db[tid]
        except gffutils.exceptions.FeatureNotFoundError:
            continue
        strand = transcript.strand
        threeputr_features = find_threeputr_features(transcript)

        for index, start, end in pool.imap_unordered(iter_mres, ((index, row['UTR_start'], row['UTR end'], strand, threeputr_features) for index, row  in group.iterrows())):
            mre_starts[index] = start - 1
            mre_ends[index] = end
            # df.loc[index, 'mre_start'] = start
            # df.loc[index, 'mre_end'] = end

assert len(mre_starts) > 0, 'db needs to find at least some of the transcript ids'
df['mre_start'] = pd.Series(mre_starts)
df['mre_start'] = df['mre_start'].fillna(-1).astype(int)
df['mre_end'] = pd.Series(mre_ends)
df['mre_end'] = df['mre_end'].fillna(-1).astype(int)

df.to_csv(snakemake.output[0])
