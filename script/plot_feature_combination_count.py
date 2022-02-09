'''
Number of identified targets when integrating different kinds of features.

Only those genes that are expressed in WT or any of the RNAi mutants with at least 1 TPM are being considered.

TargetScan designed its context++ score only for 3' UTRs of mRNAs. To prevent the exclusion of potential CDS- or 5' UTR-interactions, we decided not to use the TargetScan scores for filtering interactions.
'''

import itertools
import pdb
# already filterted for minimum of 1 TPM in any mutant
from collections import defaultdict
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import upsetplot

try:
    num_genes = pd.read_csv(snakemake.input[0], index_col=list(range(4)))['0'].sort_values(ascending=False)

    if snakemake.wildcards['direction'] == '_vertical':
        orientation = 'vertical'
    else:
        orientation = 'horizontal'
    axes = upsetplot.plot(num_genes, show_counts=True, sort_by='cardinality', orientation=orientation)
    ax = axes['intersections']
    if snakemake.wildcards['direction'] == '_vertical':
        ax.set_xlabel('Detected miRNA targets')
        ax.set_xscale('log', basex=10)
    else:
        ax.set_ylabel('Detected miRNA targets')
        ax.set_yscale('log', basey=10)
    plt.tight_layout()

    plt.savefig(snakemake.output[0])
except:
    pdb.post_mortem()
    raise
