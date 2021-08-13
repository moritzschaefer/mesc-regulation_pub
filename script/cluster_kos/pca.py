# also used for siPOOL, so should probably be moved...
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.data import ensembl_release
from moritzsphd.plot import pcaplot

df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0)

sorted_colors = [snakemake.params['sample_colors'][m[:m.rfind('_')].rstrip('KO')]
                 for m in df.columns]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    'Sample colors', sorted_colors,
    len(sorted_colors))

if snakemake.wildcards['subset'] == 'protein_coding':
    # Now only protein_coding genes

    def _biotype(i):
        try:
            return ensembl_release.gene_by_id(i).biotype
        except:
            return None

    df = df[df.index.map(_biotype) == 'protein_coding']
elif snakemake.wildcards['subset'] != 'all':
    raise ValueError('value for subset-wildcard is not supported')
df *= (1e6/df.sum())  # CPM normalization

fig, ax = plt.subplots(figsize=(6, 3.5))
ax, explained_variance = pcaplot(df, scaling='standard', ax=ax, cmap=cmap)

# formatting
ax.get_legend().remove()
sns.despine()
ax.set_title('')
ax.set_xlabel(f'{ax.get_xlabel()} ({round(explained_variance[0] * 100)}%)')
ax.set_ylabel(f'{ax.get_ylabel()} ({round(explained_variance[1] * 100)}%)')
ax.set_xticklabels([])
ax.set_yticklabels([])

fig.savefig(snakemake.output['plot'])
