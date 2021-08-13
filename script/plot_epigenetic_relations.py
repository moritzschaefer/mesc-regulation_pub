import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2, venn3

try:
    h3k27me3_log2fc = pd.read_csv(snakemake.input['h3k27me3_log2fc'], index_col=0)
    indirect_targets = pd.read_csv(snakemake.input['indirect_targets'], squeeze=True, index_col=0)
    direct_targets = pd.read_csv(snakemake.input['direct_targets'])
except NameError:
    h3k27me3_log2fc = pd.read_csv('output/h3k27me3_log2fc.csv', index_col=0)
    indirect_targets = pd.read_csv('output/indirect_targets_protein_coding.csv', squeeze=True, index_col=0)
    direct_targets = pd.read_csv('output/mirnas/interaction_ranking_protein_coding.csv')


# reduce h3k27me3_log2fc set to protein_coding

# there are two potential reasons for a missregulation: loss and increase of H3K27me3. The log2FC of H3K27me3 sohuld be inverse proportional
# We want to identify H3K27me3 marks that are consistent across mutants. Furthermore, I'll correlate these marks with my indirect_targets and look at the proportion of targets that fit (e.g. H3K27me3 is up -> target is down). This proportion should be larger than when I compare to the direct targets (those are already explained..). Also, for direct_targets, there should just be fewer consistent H3K27me3 hits

# Drosha is anyways downregulated, (whatever that means for itself..)
h3k27me3_down = h3k27me3_log2fc.index[(h3k27me3_log2fc < -0.3).sum(axis=1) >= 3]  # the cut-offs are so arbitrary. I need some kind of negative control here!
h3k27me3_up = h3k27me3_log2fc.index[((h3k27me3_log2fc > 0.1).sum(axis=1) >= 3) & ((h3k27me3_log2fc > 0.3).sum(axis=1) >= 2)]  # because of Drosha, a bit more complicated..

# there is mostly a consistent loss of H3K27me3 marks! This corresponds to what we see in Western. We hypothesize that the loss of AGO2 in the mutants contributes to the loss of H3K27me3..

# There is no overlap between H3K-ups and gene-ups. This is expected (they should be inverse proportional)
# venn3((set(h3k27me3_up), set(indirect_targets.index[indirect_targets['direction'] == 'up']), set(direct_targets['Geneid'])), ('h3k27me3_up', 'indirect_up', 'direct(_up)'))

# When comparing H3K-downs with gene-ups we do see an overlap with up-regulated genes. There is no difference between direct and indirect miRNA-targets. There is no reason why direct miRNA targets shouldn't be affected differentially methylated, so this is fine.
# A great proportion of the indirect-upregulated genes seems to be explained by the loss of H3K27me3

venn3((set(h3k27me3_down), set(indirect_targets.index[indirect_targets['direction'] == 'up']) | set(direct_targets['Geneid']), set(direct_targets['Geneid'])), ('Genes with consistently reduced H3K27me3 levels', 'Consistently upregulated genes', 'Direct miRNA targets'))

plt.savefig(snakemake.output['fig_h3k27me3_down'])

# venn2((set(h3k27me3_down), set(indirect_targets.index[indirect_targets['direction'] == 'down'])), ('h3k27me3_down', 'indirect_down'))

# On the other hand, at some genes, we also see an increase of H3K27me3, however, these only explain very few genes.
plt.figure()
venn2((set(h3k27me3_up), set(indirect_targets.index[indirect_targets['direction'] == 'down'])), ('Genes with consistently increased H3K27me3 levels', 'Consistently downregulated targets/indirect miRNA down-targets'))


plt.savefig(snakemake.output['fig_h3k27me3_up'])
