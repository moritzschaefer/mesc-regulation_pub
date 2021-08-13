import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0)
df = df[(df > 1).any(axis=1)]

corr = df.corr()

fig, ax = plt.subplots()

cmap = sns.diverging_palette(230, 20, as_cmap=True)

sns.heatmap(corr, cmap=cmap, square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax=ax)

fig.savefig(snakemake.output[0])
