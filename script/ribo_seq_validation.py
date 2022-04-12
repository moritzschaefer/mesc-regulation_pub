import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind

log = open(snakemake.log[0], 'w')

df = pd.read_excel(snakemake.input['ribo_seq'], index_col=[0, 1], skiprows=2)
df.rename(columns={'AGO2&1_log2FC': 'AGO12_log2FC', 'AGO2&1_padj': 'AGO12_padj'}, inplace=True)
mutants = snakemake.params['mutants']
ribo_log2fc = df[[m.upper() + '_log2FC' for m in mutants]]
ribo_log2fc.columns = mutants
ribo_log2fc.reset_index(level=1, inplace=True, drop=True)

control = ribo_log2fc.stack().dropna()
control = control[control.apply(type) == float]

mirna_targets = pd.read_csv(snakemake.input['interaction_ranking'], index_col=0)
ribo_selection = ribo_log2fc.reindex(mirna_targets.index)
intersection_size = len(ribo_log2fc.index.intersection(mirna_targets.index))
log.write(f'{intersection_size} of {len(mirna_targets)} miRNA targets have been observed in Ribo-seq ({round(100.0 * intersection_size / len(mirna_targets))}%)\n')

fig, ax = plt.subplots(figsize=(5, 3.3))
# ax.grid(True)
ax.axhline(y=0.5, color='gray')
ax.axvline(x=0, color='gray')
# print(ribo_selection['Dgcr8'].isna().value_counts()) # 645 vs 31
sns.ecdfplot(control.values, ax=ax, label='control (all genes)', color='black')
for mutant in mutants:
    sns.ecdfplot(ribo_selection[mutant].dropna(), ax=ax, label=mutant, color=snakemake.params['sample_colors'][mutant])
ax.set_title(f'Ribo-seq log2FC of predicted miRNA targets ({len(ribo_selection[mutant].dropna())} genes)')
ax.set_xlim([-2, 3.5])
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower right', ncol=1)
sns.despine()
ax.set_xlabel('log2FC')
# ax.grid(True)
fig.savefig(snakemake.output['plot'])

ribo_dropped = ribo_selection.dropna()
for i in range(1, (len(mutants) + 1)):
    stat = ((ribo_dropped > 0).sum(axis=1) >= i).value_counts()
    log.write(
        f'For {stat[True]} of {len(ribo_dropped)} genes ({round(100.0 * stat[True]/len(ribo_dropped))}%) we see a positive log2FC in at least {i} mutants.\n')

for mutant in mutants:
    # t-test of miRNA-targets vs control (all genes)
    control = ribo_log2fc[mutant].dropna()
    control = control[control.apply(type) == float]  # filter a couple of '-' entries -.-
    log.write(f't-test for {mutant}: {ttest_ind(ribo_selection[mutant].dropna(), control)}\n')

    log.write(f'Number of positives in {mutant}: {(ribo_selection[mutant].dropna() > 0).sum()}/{len(ribo_selection[mutant].dropna())}: {(ribo_selection[mutant].dropna() > 0).sum()/len(ribo_selection[mutant].dropna())}\n')

log.close()
