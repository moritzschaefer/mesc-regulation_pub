import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

df = pd.concat(
    [pd.DataFrame({'source': tf, 'targets': pd.read_csv(tf_file, index_col=0)['gene_id'].tolist()}) for tf, tf_file in
     zip(snakemake.params.tfs, snakemake.input.targets)])

# filter for TFs and convert to gene symbols
tf_df = pd.read_csv(snakemake.input.mirna_targets, index_col=0).query('is_tf')
df = df.loc[df['targets'].isin(tf_df.index)]
df['targets'] = tf_df.loc[df['targets']]['Gene name'].values

# substitute target gene names with our TF samples. If there is a replicate, duplicate the target
new_df = df.copy()
for tf_long in snakemake.params.tfs:
    if '_' in tf_long:
        tf = tf_long.split('_')[0]
        new_df.drop(new_df.index[new_df['targets'] == tf], inplace=True)

        tmp = df.loc[df['targets'] == tf, 'targets'].copy()
        tmp['targets'] = tf_long
        new_df = new_df.append(tmp)

# plot
graph = nx.convert_matrix.from_pandas_edgelist(df, 'source', 'targets', edge_attr=None, create_using=nx.DiGraph())
pos = nx.layout.circular_layout(graph)

nx.draw(graph, pos=pos, with_labels=True, font_weight='bold')

plt.savefig(snakemake.output[0])
