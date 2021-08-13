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

pcaplot(mirna_data, scaling='standard', ax=ax, cmap=cmap)

# formatting
ax.get_legend().remove()
sns.despine()
ax.set_title('')
ax.set_xticklabels([])
ax.set_yticklabels([])

fig.savefig(snakemake.output[0])
