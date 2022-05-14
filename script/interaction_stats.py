import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

interactions = pd.read_csv(snakemake.input.interaction_ranking)
number_gene_interactions_count = interactions.groupby('Geneid').count()['Gene name']
fig, ax = plt.subplots(figsize=(1.5, 2.7))
bin_thres = snakemake.params['mre_count_bins'] + [max(number_gene_interactions_count) + 1]
gene_groups = {f'{b1}-{b2-1}': number_gene_interactions_count.index[
    (number_gene_interactions_count >= b1) & (number_gene_interactions_count < b2)]
               for b1, b2 in zip(bin_thres[:-1], bin_thres[1:])}
hist = np.histogram(number_gene_interactions_count, bins=bin_thres)
ax.bar(list(range(len(hist[0]))), hist[0], color='black')
# pd.Series(index=plot_series.values, data=plot_series.index).plot()
# ax.set_title('Interaction count distribution')
ax.set_ylabel('Number of (target) genes')
ax.set_xlabel('Number of interactions')
# ax.set_xlim([0, 116])
# ax.set_ylim([0, 690])
sns.despine()
ax.set_xticks(list(range(len(hist[0]))))
ax.set_xticklabels(gene_groups.keys(), rotation=30, ha='right')
ax.grid(None)
plt.tight_layout()
fig.savefig(snakemake.output.interaction_histo)

#### Gene expression
fig, ax = plt.subplots(figsize=(2, 1.7))

tpm_expression = pd.read_csv(snakemake.input['mrna_data'], header=[0, 1], index_col=[0, 1])[('WT', 'tpm_expression')]

sns.distplot(np.log10(tpm_expression.loc[tpm_expression > 0.5]), color='black', ax=ax, label='expressed genes')
sns.distplot(np.log10(interactions.groupby('Geneid')['WT mRNA expression'].first()), color='#Ff8c69', ax=ax, label='Predicted miRNA targets')  # salmon1
fig.legend()
ax.set_xticks(np.log10([1, 10, 100, 1000, 10000]))
ax.set_xticklabels([1, 10, 100, 1000, 10000])
ax.set_xlim([-1, 3.8])
ax.set_xlabel('mRNA tpm')
ax.set_title('target expression (WT RNA-seq)')
sns.despine()
fig.savefig(snakemake.output.gene_expression)

##### mRNA sections
fig, ax = plt.subplots(figsize=(1.2, 2))

ax = interactions.groupby('gene_location').count()['Geneid'].plot(kind='bar', color='black')
ax.set_ylabel('Number of filtered interactions')
ax.grid(axis='x')
sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output.mrna_sections)

##### MRE type

fig, ax = plt.subplots(figsize=(1.2, 2))

ax = interactions.groupby('MRE type').count()['Geneid'].plot(kind='bar', color='black')
ax.set_ylabel('Number of filtered interactions')
ax.grid(axis='x')
sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output.mre_type)
