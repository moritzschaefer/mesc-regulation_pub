import pandas as pd

df = pd.read_csv(snakemake.input.targets, index_col=0)
mrna_df = pd.read_csv(snakemake.input.mrna_data,
                      index_col=[0, 1], header=[0, 1])
mrna_df.reset_index(drop=True, level=1, inplace=True)
l2fcs = mrna_df.xs('log2FoldChange', axis=1, level=1)[snakemake.params.mutants]

# log2FC > 0 in at least 2 mutants
# AND sum of log2fc across mutants > 0.5

if snakemake.wildcards.tf in snakemake.params.repressors:
    selected = df[df.apply(lambda row:
                           row.gene_id in l2fcs.index and
                           (l2fcs.loc[row.gene_id] < 0.0).sum() >= 2 and
                           l2fcs.loc[row.gene_id].sum() < -0.5,
                           axis=1
                           )]
else:
    selected = df[df.apply(lambda row:
                           row.gene_id in l2fcs.index and
                           (l2fcs.loc[row.gene_id] > 0.0).sum() >= 2 and
                           l2fcs.loc[row.gene_id].sum() > 0.5,
                           axis=1
                           )]

selected.to_csv(snakemake.output.targets, index=True)
