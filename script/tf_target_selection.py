import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kstest

from moritzsphd.data import genename_to_id
from moritzsphd.plot import colored_text

# pylint:disable=invalid-name,used-before-assignment
# Load inputs/params
mutants = snakemake.params['mutants']
tfs = snakemake.params['tfs']
distances = [str(d) for d in snakemake.params['distances']]
q_cutoffs = [str(q) for q in snakemake.params['q_cutoffs']]

# load mrna_expression
mrna_df = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1])
mrna_df = mrna_df[snakemake.params['mutants'] + ['WT']]
mrna_df.reset_index(drop=True, level=1, inplace=True)
# filter for expression
mrna_df = mrna_df[(mrna_df.xs('tpm_expression', axis=1, level=1) > 0.5).any(axis=1)]
all_l2fcs = mrna_df.xs('log2FoldChange', axis=1, level=1)
all_l2fcs = all_l2fcs[snakemake.params['mutants']]  # ignore WT, Ago1 and Ago2

fig, axes = plt.subplots(len(tfs), 5, squeeze=False, figsize=(21, len(tfs) * 3))

mean_dist_shift = np.ndarray((len(tfs), len(distances), len(q_cutoffs)), dtype=np.float)
median_dist_shift = np.ndarray((len(tfs), len(distances), len(q_cutoffs)), dtype=np.float)
max_dist_shift = np.ndarray((len(tfs), len(distances), len(q_cutoffs)), dtype=np.float)
peak_count = np.ndarray((len(tfs), len(distances), len(q_cutoffs)), dtype=np.float)
target_count = np.ndarray((len(tfs), len(distances), len(q_cutoffs)), dtype=np.float)


# TODO for now combine the 4 mutants! (is there another way?)
def _compute_distribution_shift(df, diff_fn=lambda targets, contr: kstest(targets, contr)[0]):
    targets = df['gene_id']
    # genes that were filtered because of low expressin (<0.5) are because they are not protein_coding (depending on the input gene set) are ignored.
    target_l2fcs = all_l2fcs.loc[all_l2fcs.index.intersection(targets)]

    try:
        # return target_l2fcs.median() - control_fcs.median()
        return diff_fn(target_l2fcs, all_l2fcs)  # probably bad to use kstest: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4342197/
    except ValueError:  # if #(targets) = 0
        return 0

def plot_distributions(control, target_dists, ax):
    MUTANT_COLORS = ['#008b8b', '#8a2be2', '#Ff8c00', '#Daa520']
    sns.kdeplot(control.stack(), cumulative=True, ax=ax, label='All genes', color='black', legend=False, gridsize=300)
    for mutant, color in zip(mutants, MUTANT_COLORS):
        sns.kdeplot(target_dists[mutant], cumulative=True, ax=ax, label=f'TF targets in {mutant}', color=color, legend=False)

    # ax.set_title(f'({len(target_dists)} targets)')  # TODO could be deleted..

    gene_id = genename_to_id(tf.split('_')[0])
    tf_l2fcs = [f'{mutant[:2]}: {all_l2fcs.loc[gene_id, mutant]:.2f}' for mutant in mutants]

    colored_text(0.2, 0.05, ['TF log2(diff. exp.):'] + tf_l2fcs, ['black'] + MUTANT_COLORS, size=8,
                 va="baseline", ha="left", multialignment="left", bbox=dict(fc="none"), ax=ax)

    ax.text(-0.9, 0.9, f'#(targets)={len(target_dists)}')

    ax.set_xlim([-1, 1])

# load heatmap data
for f in snakemake.input['csvs']:
    tf, distance, q_cutoff = re.search('targets/([^/]+)_d([0-9]+)_q([0-9\.e\-]+).csv', f).groups()
    tf_index = list(tfs.keys()).index(tf)
    distance = distance
    q_cutoff = q_cutoff
    df = pd.read_csv(f)

    for fn, hm in zip([lambda t, c: t.stack().mean() - c.stack().mean(),
                       lambda t, c: t.stack().median() - c.stack().median(),
                       lambda t, c: kstest(t.stack(), c.stack())[0],
                       lambda t, c: len(t)],
                      [mean_dist_shift, median_dist_shift,
                       max_dist_shift, target_count]):
        hm[
            tf_index,
            distances.index(distance),
            q_cutoffs.index(q_cutoff)
        ] = _compute_distribution_shift(df, fn)
    peak_count[
        tf_index,
        distances.index(distance),
        q_cutoffs.index(q_cutoff)
    ] = len(df.peak_name.drop_duplicates())

mean_dist_shift = np.nan_to_num(mean_dist_shift)
median_dist_shift = np.nan_to_num(median_dist_shift)

# the function was visually derived
score = np.abs(mean_dist_shift) * ((np.log(target_count) / np.log(1.01)) - 300)  # works reasonably well
score = np.nan_to_num(score)
score[score < -1] = -1

# plot!
selected_parameters = pd.DataFrame(index=tfs.keys(), columns=['q_cutoff', 'distance', 'mode'])

use_manual_selection = True

for tf_i, (tf, tf_meta) in enumerate(tfs.items()):
    if use_manual_selection:
        best_distance = str(tf_meta['distance'])
        best_q_cutoff = tf_meta['q']
        max_index = distances.index(best_distance), q_cutoffs.index(best_q_cutoff)
    else:
        # select the best parameters
        max_index = np.unravel_index(np.abs(score[tf_i, :, :]).argmax(), score[tf_i, :, :].shape)
        best_distance = distances[max_index[0]]
        best_q_cutoff = q_cutoffs[max_index[1]]

    # plot heatmaps
    for x_i, (ax, hm, title) in enumerate(zip(
                axes[tf_i, :4],
                [mean_dist_shift, peak_count, target_count, score],  # median_dist_shift, max_dist_shift, score
                ['mean dist. diff.', '#(peaks)', '#(targets)', 'Score f(#(targets), dist. diff.)'])):
        if x_i == 2:
            sns.heatmap(hm[tf_i, :, :], ax=ax, center=400, vmax=min(5000, hm[tf_i, :, :].max()))
        else:
            sns.heatmap(hm[tf_i, :, :], ax=ax)

        ax.set_xticks([x+0.8 for x in range(len(q_cutoffs))])
        ax.set_xticklabels(q_cutoffs, rotation=35, ha='right')
        ax.set_yticks([y+0.5 for y in range(len(distances))])
        ax.set_yticklabels(distances, rotation=0)
        if tf_i == 0:
            ax.set_title(title)
        if tf_i == 0:
            ax.set_xlabel('Peak cutoff (q-value)')
        if x_i == 1:
            ax.set_ylabel('Max. gene<->peak distance')
        if x_i == 0:
            ax.set_ylabel(f'TF: {tf}')

        # draw a cross in score, mean and #(targets)
        ax.scatter([max_index[1] + 0.5], [max_index[0] + 0.5],
                   color='red', marker='*', s=100)

    # plot selected distribution
    # get distribution for best parameters

    fn = f'output/tfs/targets/{tf}_d{best_distance}_q{best_q_cutoff}.csv'
    best_targets = pd.read_csv(fn)
    best_l2fcs = all_l2fcs.loc[all_l2fcs.index.intersection(best_targets['gene_id'])]

    ax_shift = axes[tf_i, 4]
    plot_distributions(all_l2fcs, best_l2fcs, ax=ax_shift)
    # if tf_i == 0:
    ax_shift.set_xlabel('log2(fold-change)')
    if tf_i == 0:
        ax_shift.set_title('log2fc of selected targets')
    ax_shift.set_ylabel('cumulative fraction')
    # ax_shift.legend(loc='upper left')  # too large...

    # store parameters and targets
    selected_parameters.loc[tf, 'distance'] = best_distance
    selected_parameters.loc[tf, 'q_cutoff'] = best_q_cutoff
    selected_parameters.loc[tf, 'mode'] = 'activation' if mean_dist_shift[tf_i, max_index[0], max_index[1]] > 0 else 'repression'

    best_targets.to_csv(snakemake.output.targets[tf_i])


handles, labels = ax_shift.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower right', ncol=3)

selected_parameters.to_csv(snakemake.output['parameters'], index=True)

plt.tight_layout()

fig.savefig(snakemake.output['plot'])
