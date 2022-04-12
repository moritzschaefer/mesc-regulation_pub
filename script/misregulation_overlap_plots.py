import itertools

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn3

df = pd.read_csv(snakemake.input.mrna_data, header=[0, 1], index_col=[0, 1])
# df = df[snakemake.params['mutants']]
# df = df.loc[(df.xs('tpm_expression', axis=1, level=1) > 1).any(axis=1)]  # filter lowly expression genes
padj = df.xs('padj', axis=1, level=1)
log2fc = df.xs('log2FoldChange', axis=1, level=1)

mutants = snakemake.params.mutants

# up
ups = ((padj < snakemake.params['combined_padj_threshold']) & (log2fc > snakemake.params['log2fc_threshold']))
up_sets = [set(ups.index[ups[mutant]]) for mutant in mutants]
c = venn3(up_sets, mutants)
plt.savefig(snakemake.output.up)
# c.get_patch_by_id('11').set_color('magenta')
# c.get_patch_by_id('11').set_edgecolor('none')
# c.get_patch_by_id('11').set_alpha(0.4)
plt.subplots()

# down
downs = ((padj < snakemake.params['combined_padj_threshold']) & (log2fc < snakemake.params['log2fc_threshold']))
down_sets = [set(downs.index[downs[mutant]]) for mutant in mutants]
c = venn3(down_sets, mutants)

plt.savefig(snakemake.output.down)
