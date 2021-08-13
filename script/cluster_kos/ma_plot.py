import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.plot import maplot

group = snakemake.wildcards['group']

palette = {'up': '#a00000', 'down': '#36648b', 'unchanged': 'darkgray'}
fig, ax = plt.subplots(figsize=(4.7, 3.2))

# expr = pd.read_csv(snakemake.input['expr'], sep='\t', index_col=0)
# expr *= (1e/expr.sum())

df = pd.read_excel(snakemake.input['data'], sheet_name='Main', skiprows=2, index_col=[0, 1], header=[0, 1])
try:
    df = df.xs(group, axis=1, level=0)
except KeyError:  # hacky
    df = df.xs(f'si_{group}', axis=1, level=0)
# if 'KO' in group:
#     df['mean'] = expr[['WT_1', 'WT_2', f'{group[:group.find("_")]}KO_1', f'{group[:group.find("_")]}KO_2']].mean(axis=1)
# else:
#     df['mean'] = expr[['WT_1', 'WT_2', f'{group}KO_1', f'{group}KO_2']].mean(axis=1)
maplot(data=df[df['baseMean'] > snakemake.params['baseMean']], mean='baseMean', p_threshold=snakemake.config['padj_threshold'], palette=palette)
# ax.set_title(f'{group}KO vs WT')
ax.set_title(f'{group}')

red_patch = mpatches.Patch(color=palette['up'], label='Significanly up-regulated genes')
blue_patch = mpatches.Patch(color=palette['down'], label='Significanly down-regulated genes')
plt.legend(handles=[red_patch, blue_patch], bbox_to_anchor=(0.2, -0.2), ncol=2)

sns.despine()
fig.savefig(snakemake.output[0])
