import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.plot import pcaplot

# rename is for compatibility with publication table
mirna_data = pd.read_excel(snakemake.input['mirna_data'], sheet_name='Normalized read counts', index_col=0, skiprows=2).rename(columns={'Ago2&1_1': 'Ago12_1', 'Ago2&1_2': 'Ago12_2'})[snakemake.params['samples']]

fig, ax = plt.subplots(figsize=(6, 3.5))

sorted_colors = [snakemake.params['sample_colors'][m[:m.find('_')]]
                 for m in snakemake.params['samples']]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    'Sample colors', sorted_colors,
    len(snakemake.params['sample_colors']))

ax, explained_variances = pcaplot(mirna_data, scaling='log2', ax=ax, cmap=cmap)

ax.set_xlabel(f'Principal Component 1 ({int(explained_variances[0] * 100)}%)')
ax.set_ylabel(f'Principal Component 2 ({int(explained_variances[1] * 100)}%)')

# formatting
ax.get_legend().remove()
sns.despine()
ax.set_title('')
# ax.set_xticklabels([])
# ax.set_yticklabels([])

fig.savefig(snakemake.output[0])
