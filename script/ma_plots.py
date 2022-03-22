import math

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.plot import maplot

try:
    mutants = snakemake.params['mutants']
    mrna_file = snakemake.input[0]
    outfile = snakemake.output[0]
except NameError:
    mutants = ['WT_2i']
    mrna_file = './output/mrna_data_protein_coding.csv'
    outfile = '/tmp/outfile.svg'

df = pd.read_csv(mrna_file, index_col=[0, 1], header=[0, 1])
y = math.ceil(len(mutants)/2)
x = min(2, len(mutants))
fig, axes = plt.subplots(y, x, figsize=(6*x, 4*y), sharex=True, sharey=True, squeeze=False)

palette = {'up': '#a00000', 'down': '#36648b', 'unchanged': 'darkgray'}
wt_tpm = df[('WT', 'tpm_expression')]

for i, mutant in enumerate(mutants):
    ax = axes[i//2, i%2]
    mutant_df = df.xs(mutant, axis=1, level=0).copy()
    # mutant_df['-log10(adjusted pvalue)'] = -np.log10(mutant_df['padj'])

    mutant_df['mean TPM'] = (wt_tpm + mutant_df['tpm_expression']) / 2

    #volcanoplot(data=df[df.baseMean > 20], ax=ax, alpha=0.3, show_n=True)

    maplot(data=mutant_df[mutant_df['mean TPM'] > 0.1], ax=ax, mean='mean TPM', p_threshold=snakemake.params['padj_threshold'], palette=palette)
    ax.set_title(f'{mutant}vsWT')
    if i == 0:
        red_patch = mpatches.Patch(color=palette['up'], label='Significanly up-regulated genes')
        blue_patch = mpatches.Patch(color=palette['down'], label='Significanly down-regulated genes')
        plt.legend(handles=[red_patch, blue_patch], bbox_to_anchor=(0.2, -0.2), ncol=2)

sns.despine()
fig.savefig(outfile)

# adjust position of n-labels
# for ax in axes.flatten():
#     x_low = ax.get_xlim()[0] * 0.9
#     x_high = ax.get_xlim()[1] * 0.9
#     y = ax.get_ylim()[1] * 0.9
#     texts = [t for t in ax.get_children() if 'n=' in getattr(t, 'get_text', lambda: '')()]
#     texts[0].set_position((x_low, y))
#     texts[1].set_position((x_high, y))
#     texts[1].set_horizontalalignment('right')
