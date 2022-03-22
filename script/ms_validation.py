import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.data import gn2id
from scipy.stats import ttest_ind

log = open(snakemake.log[0], 'w')

mutants = snakemake.params['mutants']
ms_log2fc = pd.read_excel(snakemake.input['ms'])

ms_log2fc = ms_log2fc[~ms_log2fc['symbol'].isna()]  # delete 9 proteins without symbol
ms_log2fc['gene_id'] = ms_log2fc['symbol'].map(gn2id)
ms_log2fc.set_index('gene_id', inplace=True)
ms_log2fc.drop(columns=['ID', 'symbol'], inplace=True)

ms_log2fc = ms_log2fc[[m.upper() + '_log2FC' for m in mutants]]
ms_log2fc.columns = mutants

# two genes had been found twice (one gene can build distinct proteins?). Since it's so few, we simple get the mean..
ms_log2fc = ms_log2fc.groupby(ms_log2fc.index).mean()

control = ms_log2fc.stack().dropna()
control = control[control.apply(type) == float]

mirna_targets = pd.read_csv(snakemake.input['interaction_ranking'], index_col=0)
ms_selection = ms_log2fc.reindex(mirna_targets.index)
intersection_size = len(ms_log2fc.index.intersection(mirna_targets.index))
log.write(f'{intersection_size} of {len(mirna_targets)} miRNA targets have been observed in MS ({round(100.0 * intersection_size / len(mirna_targets))}%)\n')

fig, ax = plt.subplots(figsize=(5, 3.5))
# ax.grid(True)
ax.axhline(y=0.5, color='gray')
ax.axvline(x=0, color='gray')
# print(ms_selection['Dgcr8'].isna().value_counts()) # 645 vs 31
sns.ecdfplot(control.values, ax=ax, label='control (all genes)', color='black')
for mutant in mutants:
    sns.ecdfplot(ms_selection[mutant].dropna(), ax=ax, label=mutant, color=snakemake.params['sample_colors'][mutant])
ax.set_title(f'MS log2FC of predicted miRNA targets ({len(ms_selection[mutant].dropna())} genes)')
ax.set_xlim([-1.5, 2])
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower right', ncol=1)
sns.despine()
ax.set_xlabel('log2FC')
fig.savefig(snakemake.output['plot'])

ms_dropped = ms_selection.dropna()
for i in range(1, (len(mutants) + 1)):
    stat = ((ms_dropped > 0).sum(axis=1) >= i).value_counts()
    log.write(
        f'For {stat[True]} of {len(ms_dropped)} genes ({round(100.0 * stat[True]/len(ms_dropped))}%) we see a positive log2FC in at least {i} mutants.\n')

for mutant in mutants:
    # t-test of miRNA-targets vs control (all genes)
    log.write(f't-test for {mutant}: {ttest_ind(ms_selection[mutant].dropna(), ms_log2fc[mutant].dropna())}\n')

    log.write(f'Number of positives in {mutant}: {(ms_selection[mutant].dropna() > 0).sum()}/{len(ms_selection[mutant].dropna())}: {(ms_selection[mutant].dropna() > 0).sum()/len(ms_selection[mutant].dropna())}\n')

log.close()
