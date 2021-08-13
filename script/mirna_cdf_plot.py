import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_csv(snakemake.input[0], header=[0, 1], index_col=0)
mutants = snakemake.params['mutants']

# filter for minimal expression
df = df[df[('WT', 'Expression')] > 1]
df = df.xs('log2FoldChange', axis=1, level=1)[mutants]

fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))

ax.axhline(y=0.5, color='gray')
ax.axvline(x=0, color='gray')

plot_df = df.melt(value_vars=mutants, ignore_index=False, var_name='mutant', value_name='log2(fold-change)')
sns.ecdfplot(plot_df, hue='mutant', x='log2(fold-change)', palette=snakemake.params['sample_colors'])

sns.despine()
plt.tight_layout()
fig.savefig(snakemake.output[0])
