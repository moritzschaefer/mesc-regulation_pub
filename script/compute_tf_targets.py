import gffutils
import numpy as np
import pandas as pd

db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)
BED_NAMES = [
    'chrom', 'start', 'end', 'peak_name', 'score', 'strand',
    'signalValue', 'pvalue', 'qvalue', 'peak_offset',

    'transcript_chrom', 'transcript_source', 'transcript_type',
    'transcript_start', 'transcript_end', 'dunno1', 'transcript_strand',
    'dunno2', 'attributes',
    'transcript_orf_start', 'transcript_orf_end', 'dunno3']

df = pd.read_csv(snakemake.input.bed, header=None, sep='\t', names=BED_NAMES)
# find the distance between peak and gene
# if d1 and d2 are negative, then distance = 0
# else, the positive one is the closer one

def _distance(row):
    d1 = row['start'] - row['transcript_end']
    d2 = row['transcript_start'] - row['end']
    if d1 <= 0 and d2 <= 0:
        return 0
    else:
        return max(d1, d2)

# we define gene-start location as TSS
def _gene_start_distance(row):
    start = row['transcript_start'] if row['transcript_strand'] == '+' else row['transcript_end']
    return min(abs(start - row['start']), abs(start - row['end']))

df['distance'] = df.apply(_distance, axis=1)
df['tss_distance'] = df.apply(_gene_start_distance, axis=1)

q_cutoff = float(snakemake.wildcards['q_cutoff'])

# filter for gene-features with sufficiently high cutoffs (peak) and max distance
df = df.loc[(df.transcript_type == 'gene') &
            (df['qvalue'] > -np.log10(q_cutoff)) &
            (df['distance'] < int(snakemake.wildcards['distance']))]

# get gene_id and reduce to one entry per gene
attributes = pd.DataFrame(df['attributes'].apply(
    lambda v: dict([s.strip('"') for s in a.split(' ')]
                   for a in v.rstrip(';').split('; '))).tolist(),
                          index=df.index)
try:
    df['gene_id'] = attributes['gene_id']
except KeyError:
    df['gene_id'] = []

df = df.groupby('gene_id').first().reset_index()  # we only need each gene once. 

# peak-associated values are a bit non-sense, since I am taking the "first" peak for each gene and ignore the others
df[['chrom', 'start', 'end', 'peak_name', 'score', 'strand', 'signalValue',
    'pvalue', 'qvalue', 'peak_offset', 'gene_id', 'distance', 'tss_distance']]\
    .to_csv(snakemake.output[0], index=False)
