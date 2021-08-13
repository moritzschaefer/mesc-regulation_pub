import matplotlib.pyplot as plt
import pandas as pd
import upsetplot

data = upsetplot.from_contents(
    {tf: targets for tf, targets in
     zip(snakemake.params.tfs,
         [pd.read_csv(tf, index_col=0)['gene_id'].tolist() for tf in snakemake.input.targets])
     })

# only take the top 15 overlaps
selection = pd.DataFrame(data=data.index)[0].value_counts()[:25]
upsetplot.plot(data.loc[data.index.isin(selection.index)],
               sort_by='cardinality', show_counts='%d')

plt.savefig(snakemake.output[0])
