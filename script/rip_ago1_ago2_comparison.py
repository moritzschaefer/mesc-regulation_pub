import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_excel(snakemake.input['mirna_data'], skiprows=2, index_col=0)

# df = np.log2(df+1)

fig, ax = plt.subplots(figsize=(3.5, 3.5))

ax.plot([0, df[['RIP_AGO2', 'RIP_AGO1']].max().max()], [0, df[['RIP_AGO2', 'RIP_AGO1']].max().max()], color='gray')
sns.scatterplot(data=df, x='RIP_AGO1', y='RIP_AGO2', color='black')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_yticks([10, 1000, 100000])
sns.despine()

fig.savefig(snakemake.output[0])
