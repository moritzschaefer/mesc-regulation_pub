import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3

mtb_df = pd.read_excel(snakemake.input['mtb_mmu'][0])

mirna_loading = pd.read_excel(snakemake.input['mirna_data'], skiprows=2, index_col=0)[['RIP_AGO2', 'RIP_AGO1']].mean(axis=1)
mirnas = mirna_loading.index[mirna_loading > snakemake.params['mirna_threshold']]
target_genes_all = mtb_df.loc[mtb_df['miRNA'].isin(mirnas), 'Target Gene'].drop_duplicates()

target_genes_confident = mtb_df.loc[mtb_df['miRNA'].isin(mirnas) & mtb_df['Experiments'].str.contains('//'), 'Target Gene'].drop_duplicates()


interactions = pd.read_csv(snakemake.input['interaction_ranking'])

venn3([set(l.values) for l in [target_genes_all, target_genes_confident, interactions['Gene name'].drop_duplicates()]], ('miRTarBase full', 'miRTarBase confidence', 'Integrative analysis'))

plt.savefig(snakemake.output[0])
