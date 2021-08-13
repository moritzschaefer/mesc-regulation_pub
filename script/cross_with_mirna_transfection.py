import numpy as np
import pandas as pd

# The studied miRNAs are (Liu et al. 2019):
old_names = 'hsa-miR-126-3p hsa-let-7c-5p hsa-miR-155-5p hsa-miR-206 hsa-miR-126-5p hsa-miR-107 hsa-miR-15a-5p hsa-miR-210-3p hsa-miR-146a-5p hsa-miR-10a-5p hsa-miR-16-5p hsa-miR-21-5p hsa-miR-200a-3p hsa-miR-124-3p hsa-miR-17-5p hsa-miR-31-5p hsa-miR-9-3p hsa-miR-133b hsa-miR-200b-3p hsa-miR-34a-5p hsa-miR-9-5p hsa-miR-142-3p hsa-miR-200c-3p NC hsa-miR-145-5p hsa-miR-193b-3p'.split(' ')
new_names = 'mmu-miR-126a-3p mmu-let-7c-5p mmu-miR-155-5p mmu-miR-206-3p mmu-miR-126a-5p mmu-miR-107-3p mmu-miR-15a-5p mmu-miR-210-3p mmu-miR-146a-5p mmu-miR-10a-5p mmu-miR-16-5p mmu-miR-21a-5p mmu-miR-200a-3p mmu-miR-124-3p mmu-miR-17-5p mmu-miR-31-5p mmu-miR-9-3p mmu-miR-133b-3p mmu-miR-200b-3p mmu-miR-34a-5p mmu-miR-9-5p mmu-miR-142a-3p mmu-miR-200c-3p NC mmu-miR-145a-5p mmu-miR-193b-3p'.split(' ') # already translated manually

name_translation = dict([(old, new) for old, new in zip(old_names, new_names)])

# for now just sum up the two replicates. We can maybe later do a proper DEG analysis
dfs = [pd.read_excel(snakemake.input['liu_data'][0], skiprows=[0, 1, 3], sheet_name=f'RNA-Seq Replicate {repl}', index_col=[26, 27, 28]) for repl in [1, 2]]
df = dfs[0] + dfs[1]

df.rename(columns=name_translation, inplace=True)

df = 1e6 * df / df.sum(axis=0)  # CPM normalization

# downregulated genes
positives = (df < 0.6 * df.median(axis=1).values.reshape(-1, 1))  # 40% reduction as described by Liu et al.

# store the targets
out_dfs = []
for col in positives:
    targets = df.loc[positives[col].values]
    sub_df = pd.DataFrame.from_dict({
        'miRNA': col,
        'target': targets.index.get_level_values('Gene Symbol'),
        'log2FoldChange': np.log2(targets[col] / targets.median(axis=1)).values,
        'control_expression': targets['NC'].values
    })

    out_dfs.append(sub_df.sort_values('log2FoldChange'))

pd.concat(out_dfs).to_csv(snakemake.output['liu_targets'])

mirna_interactions = pd.read_csv(snakemake.input['mirna_interactions'])

new_names.remove('NC')
with open(snakemake.log[0], 'w') as log:
    for mirna in new_names:
        pipeline_targets = mirna_interactions.loc[mirna_interactions.miRNA == mirna, 'Gene name'].str.upper()  # upper is necessary to be compatible with liu's annotation

        # compare ours to all others to get control score
        other_targets = []
        for mirna2 in new_names:
            # TODO we should only check genes that are expressed in my context...
            liu_mirna_targets = positives.index[positives[mirna2]].get_level_values('Gene Symbol')
            if mirna2 == mirna:
                targets = (set(liu_mirna_targets) & set(pipeline_targets))
            else:
                other_targets.extend(set(liu_mirna_targets) & set(pipeline_targets))

        log.write(f'{mirna}: {", ".join(targets)}, (other miRNAs targets, i.e. negative control targets: {", ".join(other_targets)})\n')
