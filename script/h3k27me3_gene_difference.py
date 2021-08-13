import gffutils
import numpy as np
import pandas as pd

from moritzsphd.util.bigwig import tss_means

look_region = 'tss-3000+3000'  # even more than 5000 is possible (according to madlen)

try:
    db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)
except NameError:
    db = gffutils.FeatureDB('./ref/gencode.db', keep_order=True)


if snakemake.wildcards['subset'] == 'all':
    genes = [gene for gene in db.features_of_type('gene')]
if snakemake.wildcards['subset'] == 'protein_coding':
    genes = [gene for gene in db.features_of_type('gene') if gene.attributes['gene_biotype'][0] == 'protein_coding']

h3k27me3_df = pd.DataFrame({sample: tss_means(fn, genes, db, look_region) for fn, sample in zip(snakemake.input['h3kdata'], snakemake.params['samples'])})

# h3k27me3_df = h3k27me3_df.groupby(axis=1, level=0).mean()
# take the first replicate, because the second one sux:
h3k27me3_df = h3k27me3_df.xs(1, axis=1, level=1)

h3k27me3_df.to_csv(snakemake.output['per_gene_tss'])

# the bigwigs contain log2FCs. convert log2 ratio to ratio
h3k27me3_df = 2**h3k27me3_df

wt = h3k27me3_df.pop('WT')

h3k27me3_df = h3k27me3_df.divide(wt.to_list(), axis=0)  # to_list is a performance-workaround
h3k27me3_df = np.log2(h3k27me3_df)

# inf_to_max = h3k27me3_df[~np.isinf(h3k27me3_df)].max().max()
# inf_to_min = h3k27me3_df[~np.isinf(h3k27me3_df)].min().min()
# h3k27me3_df.where(h3k27me3_df != np.inf, inf_to_max, inplace=True)
# h3k27me3_df.where(h3k27me3_df != -np.inf, inf_to_min, inplace=True)
h3k27me3_df.dropna(inplace=True)

h3k27me3_df.to_csv(snakemake.output['per_gene_tss_log2fc'])
