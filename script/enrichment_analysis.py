import matplotlib.pyplot as plt
import pandas as pd
from moritzsphd.analysis.enrich import enrichr
from moritzsphd.data import gid2n

gene_sets = ['WikiPathways_2019_Mouse', 'GO_Biological_Process_2021',
             'GO_Molecular_Function_2021', 'GO_Cellular_Component_2021',
             'BioPlanet_2019', 'KEGG_2019_Mouse']
column = snakemake.params['column']
if column:
    direct_mirna_interactions = pd.read_csv(snakemake.input[0])
    gene_list = direct_mirna_interactions[column].drop_duplicates().to_list()
elif column is None:
    with open(snakemake.input[0], 'r') as f:
        gene_list = [gid2n(gid) for gid in f.read().strip('\n').split(',')]
df, ax = enrichr(gene_list, gene_sets=gene_sets, filter_padj=snakemake.params['filter_padj'], adjust_texts=True)
df.to_csv(snakemake.output['enrichment_table'])
plt.savefig(snakemake.output['enrichment_plot'])
