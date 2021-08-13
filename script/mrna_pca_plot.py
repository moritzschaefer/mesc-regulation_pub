import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from moritzsphd.plot import pcaplot

counts = pd.read_csv(snakemake.input['expression'], sep='\t', index_col=0)
counts *= (1e6/counts.sum())
fig, ax = plt.subplots(figsize=(6, 3.5))

sorted_colors = [snakemake.params['sample_colors'][m[:m.find('_')]]
                 for m in snakemake.params['samples']]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    'Sample colors', sorted_colors,
    len(snakemake.params['sample_colors']))
pcaplot(counts, scaling='standard', ax=ax, cmap=cmap)

# formatting
ax.get_legend().remove()
sns.despine()
ax.set_title('')
ax.set_xticklabels([])
ax.set_yticklabels([])

fig.savefig(snakemake.output[0])
