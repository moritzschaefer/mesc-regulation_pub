import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

try:
    mutants = snakemake.params['mutants']
    mrna_file = snakemake.input[0]
    outfile = snakemake.output[0]
except NameError:
    mutants = ['Dgcr8', 'Drosha', 'Dicer', 'Ago12']
    mrna_file = './output/mrna_data_protein_coding.csv'
    outfile = '/tmp/outfile.png'

filtered = snakemake.wildcards['filtered'] == '_filtered'

df = pd.read_csv(mrna_file, index_col=[0, 1], header=[0, 1])
expression = df.xs('tpm_expression', axis=1, level=1)
expression = np.log2(expression + 1)
log2fc = df.xs('log2FoldChange', axis=1, level=1)
padj = df.xs('padj', axis=1, level=1)

if filtered:
    plotdf = log2fc.loc[(padj < 0.1).any(axis=1)][mutants].dropna()
else:
    plotdf = log2fc[mutants].dropna()

grid = sns.pairplot(plotdf, plot_kws=dict(marker="+", linewidth=0.5, alpha=0.2, rasterized=True))
grid.fig.set_figheight(10)
grid.fig.set_figwidth(10)

for x, y in np.ndindex(len(mutants), len(mutants)):
    ax = grid.axes[y, x]
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 10])
    mutant_x = mutants[x]
    mutant_y = mutants[y]

    pr, p = pearsonr(plotdf[mutant_x], plotdf[mutant_y])
    ax.text(-8, 8, f'r={pr:.2}')

grid.savefig(outfile)
