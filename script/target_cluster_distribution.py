import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from moritzsphd.util.mirna.clustering import cluster_to_name, mirbase_clusters


def process_count_table(count_table, clusters, show_n=10):
    '''
    Args:
        count_table: DataFrame with index name and one column count
        clusters: mirna clusters with same naming scheme as df
        show_n: Show the first n clusters
    '''

    count_table['cluster'] = ''
    for i, cluster in enumerate(clusters):
        cluster_name = cluster_to_name(cluster)
        count_table.loc[count_table.index.isin(cluster), ('cluster', )] = \
            cluster_name

    grouped_counts = count_table.groupby('cluster').sum()
    sorted_counts = grouped_counts.sort_values('count', ascending=False)

    other_clusters = sorted_counts.iloc[show_n:, :].sum()
    other_clusters.name = 'other'
    sorted_counts.drop(sorted_counts.iloc[show_n:].index, inplace=True)
    sorted_counts = sorted_counts.append(other_clusters)

    return sorted_counts


clusters = mirbase_clusters()

# *Compute cluster expression*
df = pd.read_csv(snakemake.input['mirna_data'], index_col=0, header=[0, 1])[[('WT', 'Expression')]]
df.columns = ['count']

num_clusters = 5
df = process_count_table(df, clusters, num_clusters)
df.rename(columns={'count': 'expression'}, inplace=True)

mirna_to_clusters = {m: cluster_to_name(c) for c in clusters for m in c}

# *Add Target count*
target_list = pd.read_csv(snakemake.input['interaction_ranking'])
target_list['cluster'] = target_list['miRNA'].apply(mirna_to_clusters.get)
cluster_target_counts = target_list.groupby('cluster')['Geneid'].nunique()
# compute "other"
cluster_target_counts['other'] = cluster_target_counts.sum() - cluster_target_counts.reindex(df.index).sum()
df['Interaction count'] = cluster_target_counts.loc[df.index]

# ** Pie plots **
fig, axes = plt.subplots(1, 2, figsize=(10, 10))
cluster_color_order = []
# gist_earth ar terrain
cmap = sns.color_palette('flare', num_clusters+3, as_cmap=True)

colors = []

for cluster, rest in df.iterrows():
    if cluster not in cluster_color_order:
        cluster_color_order.append(cluster)

    colors.append(cmap(cluster_color_order.index(cluster)))
ax = df.plot(
    kind='pie',
    ax=axes[0],
    y='expression',
    legend=False,
    autopct=lambda pct: str(round(pct/100.*np.sum(df['expression']))),
    startangle=90,
    labels=df.index[:num_clusters].tolist() + [''] * (len(df)-(num_clusters + 1)) + df.index[-1:].tolist(),
    colors=colors)
ax.set_ylabel('expression')

ax = df.plot(
    kind='pie',
    ax=axes[1],
    y='Interaction count',
    legend=False,
    labels=df.index[:num_clusters].tolist() + [''] * (len(df)-(num_clusters + 1)) + df.index[-1:].tolist(),
    autopct=lambda pct: str(round(pct/100.*np.sum(df['Interaction count']))),
    startangle=90,
    colors=colors)

plt.savefig(snakemake.output['pieplot'])

# ** bar plot **

plotdf = df.iloc[:-1]  # drop "other" cluster
plotdf['Target count'] = target_list.groupby('cluster')['Geneid'].nunique()
width = 0.27
fig, ax1 = plt.subplots(figsize=(5, 5))
ax2 = ax1.twinx()
colors = ['#696969', '#6ca6cd']
# colors = ['#696969', '#8b0000', '#ffffff']
columns = ['expression', 'Target count']
# columns = ['expression', 'AGOs RIP', 'Target count']
axes = [ax1, ax2]

for i, (color, column, ax) in enumerate(zip(colors, columns, axes)):
    plotdf[column].plot(kind='bar', color=color, ax=ax, width=width, position=1-i)
    ax.grid(False)

ax1.set_xlim(-0.4, len(plotdf) - 0.3) # #(clusters) - 0.3
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=30, ha='right')
ax1.set_title('Expression, #(gene targets) of top expressed clusters')
ax1.set_ylabel('expression (CPM-like)')
ax2.set_ylabel('Target count')

patches = [plt.plot([],[], color=color,
                     label="{:s}".format(text) )[0]  for color, text in zip(colors, columns)]
plt.legend(handles=patches, loc='upper right', ncol=1, numpoints=1)
plt.tight_layout()
fig.savefig(snakemake.output['barplot'])
