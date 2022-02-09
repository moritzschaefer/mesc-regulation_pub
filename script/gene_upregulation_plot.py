'''
This script was hackyly tuned to be able to also plot down genes
'''
import itertools

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_csv(snakemake.input['mrna_data'], header=[0, 1], index_col=[0, 1])
df = df[snakemake.params['mutants']]
# df = df.loc[(df.xs('tpm_expression', axis=1, level=1) > 1).any(axis=1)]  # filter lowly expression genes
padj = df.xs('padj', axis=1, level=1)
log2fc = df.xs('log2FoldChange', axis=1, level=1)

direction = snakemake.wildcards['direction']

max_padjs = [0.2, 0.1, 0.05, 0.01]
min_ups = [2, 3, 4]
plot_df = pd.DataFrame(index=max_padjs, columns=min_ups, dtype=int)
for max_padj, min_up in itertools.product(max_padjs, min_ups):
    if direction == 'up':
        count = (((padj < max_padj) & (log2fc > 0.0)).sum(axis=1) >= min_up).sum()
    else:
        count = (((padj < max_padj) & (log2fc < 0.0)).sum(axis=1) >= min_up).sum()
    plot_df.loc[max_padj, min_up] = int(count)

plot_df.columns = [str(v) for v in plot_df.columns]
plot_df.index = [str(v) for v in plot_df.index]
fig, ax = plt.subplots(figsize=(4.2, 3.5))
cmap = sns.color_palette("flare", as_cmap=True)
sns.heatmap(plot_df, annot=True, fmt='.0f', cmap=cmap, ax=ax, vmin=0)
ax.set_xlabel(f'{direction.title()}regulated in at least <x> mutants')
ax.set_ylabel('Maximum FDR')
plt.tight_layout()
plt.savefig(snakemake.output[0])
