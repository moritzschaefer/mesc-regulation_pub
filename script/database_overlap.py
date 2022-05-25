import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2

mtb_df = pd.read_excel(snakemake.input['mtb_mmu'][0])

mirna_loading = pd.read_excel(snakemake.input['mirna_data'], skiprows=2, index_col=0)[['RIP_AGO2', 'RIP_AGO1']].mean(axis=1)
mirnas = mirna_loading.index[mirna_loading > snakemake.params['mirna_threshold']]
target_genes_all = mtb_df.loc[mtb_df['miRNA'].isin(mirnas), 'Target Gene'].drop_duplicates()

interactions = pd.read_csv(snakemake.input['interaction_ranking'])
fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.7))

gene_sets = {k: set(l.values) for k, l in zip(('miRTarBase full', 'Integrative analysis'), [target_genes_all, interactions['Gene name'].drop_duplicates()])}
venn2(gene_sets.values(), gene_sets.keys(), ax=ax)

combined = list(set.union(*gene_sets.values()))
df = pd.DataFrame(index=combined,
                  data={label: [gene in gene_set for gene in combined]
                        for label, gene_set in gene_sets.items()})

df.to_excel(snakemake.output['data'])
plt.savefig(snakemake.output['plot'])
