import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3

ago2_heap_data = pd.read_csv(snakemake.input['ago2_heap_data'], index_col=0)  # of note, 6mers are filtered here!
ts_predictions = pd.read_csv(snakemake.input['targetscan_prediction'], index_col=0)


mrna_data = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1])
mrna_data = mrna_data[snakemake.params['mutants'] + ['WT']]

# get up-genes
log2fc = mrna_data.xs('log2FoldChange', axis=1, level=1)  # .loc[:, (slice(None), 'log2FoldChange')].
padj = mrna_data.xs('padj', axis=1, level=1)
expr = mrna_data.xs('tpm_expression', axis=1, level=1)
log2fc.index = log2fc.index.droplevel(1)
padj.index = padj.index.droplevel(1)
expr.index = expr.index.droplevel(1)
expressed_genes = set(expr.index[(expr > 1).any(axis=1)])

up_genes = mrna_data.index[(log2fc > snakemake.params['min_log2fc_lower']).all(axis=1) &
                           (((log2fc > snakemake.params['min_log2fc_upper']) &
                             (padj < snakemake.params['max_padj'])).sum(axis=1) >= snakemake.params['min_num_up_genes'])]


ts_genes = set(ts_predictions['Geneid'].unique()) & expressed_genes
ago2_genes = set(ago2_heap_data['gene_id'].unique()) & expressed_genes

ts_predictions = ts_predictions.loc[ts_predictions.Geneid.isin(up_genes.get_level_values(0))]

plt.subplots(figsize=(6, 6))
venn3([set(up_genes.get_level_values(0)), ts_genes, ago2_genes], ['Up-regulated genes', 'TS- predicted genes', 'AGO2-binding genes'])
plt.tight_layout()
plt.savefig(snakemake.output['overlap_plot'])
