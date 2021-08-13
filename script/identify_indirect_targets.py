import pandas as pd

mrna_data = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1])
mrna_data = mrna_data[snakemake.params['mutants'] + ['WT']]
mirna_interactions = pd.read_csv(snakemake.input['mirna_interactions'])

# identify commonly regulated genes:

log2fc = mrna_data.xs('log2FoldChange', axis=1, level=1)   #.loc[:, (slice(None), 'log2FoldChange')].
padj = mrna_data.xs('padj', axis=1, level=1)
log2fc.index = log2fc.index.droplevel(1)
padj.index = padj.index.droplevel(1)

up_genes = mrna_data.index[(log2fc > -0.2).all(axis=1) & (((log2fc > 0.5) & (padj < 0.1)).sum(axis=1) >= 3)]
down_genes = mrna_data.index[(log2fc < 0.2).all(axis=1) & (((log2fc < -0.5) & (padj < 0.1)).sum(axis=1) >= 3)]

# now, delete the up_genes that have been explained by mirnas (for now I just delete all of them. I could also filter for interaction_score first...)
index = up_genes.difference(mirna_interactions.set_index(['Geneid', 'Gene name']).index).append(down_genes)

df = pd.DataFrame(index=index)
df['direction'] = 'up'
df.loc[down_genes, 'direction'] = 'down'

# TODO maybe add co-regulation score.. (padj-combination as in networking..)

df.to_csv(snakemake.output['indirect_targets'])
