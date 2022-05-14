import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2

mtb_df = pd.read_excel(snakemake.input['mtb_mmu'][0])

mirna_loading = pd.read_excel(snakemake.input['mirna_data'], skiprows=2, index_col=0)[['RIP_AGO2', 'RIP_AGO1']].mean(axis=1)
mirnas = mirna_loading.index[mirna_loading > snakemake.params['mirna_threshold']]
target_genes_all = mtb_df.loc[mtb_df['miRNA'].isin(mirnas), 'Target Gene'].drop_duplicates()

interactions = pd.read_csv(snakemake.input['interaction_ranking'])
fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.7))

venn2([set(l.values) for l in [target_genes_all, interactions['Gene name'].drop_duplicates()]], ('miRTarBase full', 'Integrative analysis'), ax=ax)

plt.savefig(snakemake.output[0])
