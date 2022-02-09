import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

# TODO maybe add (mir)tarbase somehow

interactions = pd.read_csv(snakemake.input['unfiltered_interactions'])

# TODO normalization to [0, 1] should be somehow tested/controlled. Should I really take the log here? I think I took a look at this and it looked like a good idea
interactions['ago2_heap_binding_scaled'] = MinMaxScaler().fit_transform(
    np.log2(interactions['AGO2 HEAP peak'].fillna(0) + 1).values.reshape(-1, 1))

# TODO fillna(0) should be replaced with something binding-site dependent (e.g. better score for 8mers). Idea would be the CNN model from Bartel :cite:`briskin20_bioch_basis_cooper_action_micror`
interactions['ts_score_scaled'] = MinMaxScaler().fit_transform(
    -interactions['weighted context++ score'].fillna(0).values.reshape(-1, 1))

# negative interaction (because miRNAs..)
interactions['interaction_potential'] = (interactions['ago2_heap_binding_scaled'] + interactions['ts_score_scaled']) \
    / -2.0

# type of interaction
interactions['interaction_type'] = 'mirna'

# TODO, maybe use two separate scores for an interaction, instead of a combined one..
edges = interactions.set_index(['miRNA', 'Geneid'])[['interaction_potential', 'interaction_type']]

# nodes.to_csv(snakemake.output['nodes'])
edges.to_csv(snakemake.output['edges'])
