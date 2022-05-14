import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.plot import pcaplot

counts = pd.read_csv(snakemake.input['expression'], sep='\t', index_col=0)
counts *= (1e6/counts.sum())
fig, ax = plt.subplots(figsize=(3.5, 2.0))

sorted_colors = [snakemake.params['sample_colors'][m[:m.find('_')]]
                 for m in snakemake.params['samples']]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    'Sample colors', sorted_colors,
    len(snakemake.params['sample_colors']))
ax, explained_variances = pcaplot(counts, scaling='log2', ax=ax, cmap=cmap)

ax.set_xlabel(f'Principal Component 1 ({int(explained_variances[0] * 100)}%)')
ax.set_ylabel(f'Principal Component 2 ({int(explained_variances[1] * 100)}%)')

# formatting
ax.get_legend().remove()
sns.despine()
ax.set_title('')
# ax.set_xticklabels([])
# ax.set_yticklabels([])

fig.savefig(snakemake.output[0])
