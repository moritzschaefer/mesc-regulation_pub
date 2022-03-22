import matplotlib.pyplot as plt
import pandas as pd
from moritzsphd.analysis.enrich import enrichr

direct_mirna_interactions = pd.read_csv(snakemake.input[0])
gene_list = direct_mirna_interactions['Gene name'].drop_duplicates()
df, ax = enrichr(gene_list.tolist(), filter_padj=snakemake.params['filter_padj'], adjust_texts=True)
df.to_csv(snakemake.output['enrichment_table'])
plt.savefig(snakemake.output['enrichment_plot'])
