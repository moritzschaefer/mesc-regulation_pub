# filter for protein coding and do the same thing again

import gffutils
import numpy as np
import pandas as pd


def biotype_from_ensid(v):
    try:
        return db[v].attributes['gene_biotype'][0]
    except (IndexError, gffutils.exceptions.FeatureNotFoundError):
        return None

df = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1])
db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)

biotype = df.index.get_level_values(0).map(biotype_from_ensid)
biotype.index = df.index

protein_coding_df = df[biotype == 'protein_coding']

with open(snakemake.log[0], 'w') as f:
    f.write(f'{len(protein_coding_df)}/{len(df)} genes were protein_coding')

protein_coding_df.to_csv(snakemake.output[0])
