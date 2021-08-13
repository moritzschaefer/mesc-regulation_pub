'''
Mostly Figure2 plots for my paper
'''
from moritzsphd.util import remote_file
from moritzsphd.data import ensembl_release


def _target_region(wildcards, output):
    # except KeyError:
    #     gene = ensembl_release.genes_by_name(wildcards.gene)[0]
    #     region = {
    #         'chr': gene.contig,
    #         'start': gene.start,
    #         'end': gene.end
    #     }
    #     print(gene)
    region = config['target_params'][wildcards.gene]
    return f'{region["chr"]}:{region["start"]}-{region["end"]}'

rule heap_mres:
    '''
    Filter HEAP MREs for miRNA expression and save as bed file
    '''
    input:
        mirna_data='output/mirna_data.csv',
        ago2_heap_data='output/ago2_heap_mesc_mres.csv',
    params:
        min_mirna_expression=config['mirna_threshold'],
    output:
        'output/ago2_heap_mres.bed'
    conda:
        "../env/python.yaml"
    script:
        '../script/heap_mres.py'

rule expressed_mres_bed:
    '''
    Generate bed track for expressed MREs within CDS, 3'UTR and 5'UTR. For all...
    Since we only need this for plotting tracks, we can compute this for the sequences in config.yaml
    '''
    input:
        annotation='ref/gencode.db',
        config='config.yaml',  # if the config changes, this one should be recomputed
        mirna_data='output/mirna_data.csv',
        heap_peaks=remote_file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139345&format=file&file=GSE139345_mESC_peaks.csv.gz', gunzip=True)
    params:
        regions=config['target_params'],
        min_mirna_expression=config['mirna_threshold'],
    output:
        'output/track_plots/expressed_mres.bed'
    conda:
        "../env/python.yaml"
    script:
        '../script/expressed_mres_bed.py'

rule heap_peaks_to_bed:
    input: remote_file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139345&format=file&file=GSE139345_mESC_peaks.csv.gz', gunzip=True)
    output: 'output/track_plots/heap_peaks.bed'
    run:
        import pandas as pd
        df = pd.read_csv(input[0])
        df['seqnames'] = df['seqnames'].str.replace('chr', '')
        df['strand'] = '.'
        df[['seqnames', 'start', 'end', 'name', 'score', 'strand']].to_csv(output[0], sep='\t', index=False, header=False)

rule plot_target_track:
    'Plot the features of a given data point'
    'Genome track for putative miRNA targets'
    input:
        tracks='misc/tracks_manual.ini',
        # these ones are part of the ini-file and therefore required
        # genes=remote_file(config['gencode_gtf'], with_extension=True, gunzip=True),
        # heap=remote_file('cclab:/home/schamori/HEAP-Ago2/heap_ratio.bw'),  # generated from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139345
        # heap_peaks='output/track_plots/heap_peaks.bed',
        # heap_mres='output/ago2_heap_mres.bed',
        # expressed_mres='output/track_plots/expressed_mres.bed',
    output:
        'plot/target_data/track_{gene}.pdf',
    params:
        # region=lambda wildcards, output: f'{region["chr"]}:{region["start"]-500}-{region["end"]+500}'
        region=_target_region
    conda:
        "../env/python.yaml"
    shell: '''
        pyGenomeTracks --tracks {input.tracks} --region {params.region} --outFileName {output} --plotWidth 20 --trackLabelFraction '0.0'
    '''

rule plot_target_mirna_expr:
    input:
        mirna_data='output/mirna_data.csv',
        interaction_data='output/mirnas/interaction_ranking_all.csv',
        unfiltered_interaction_data='output/mirnas/unfiltered_interaction_ranking_all.csv'
    output:
        'plot/target_data/mirna_expr_{gene}.pdf'
    params:
        mirna_threshold=config['mirna_threshold'],
        plot_non_heaps=True,
    # conda:
    #     "../env/python.yaml"
    run:
        import matplotlib.pyplot as plt
        import pandas as pd
        import numpy as np
        import seaborn as sns

        interaction_df = pd.read_csv(input['interaction_data'])
        unfiltered_interaction_df = pd.read_csv(input['unfiltered_interaction_data']).query(f'`WT miRNA expression` > {params.mirna_threshold}')
        mirnas = interaction_df.loc[interaction_df['Gene name'] == wildcards['gene'],'miRNA']

        wt_expr = pd.read_csv(input['mirna_data'], header=[0, 1], index_col=0)[('WT','Expression')]
        wt_expr = wt_expr[wt_expr > params['mirna_threshold']]
        log_expr = np.log2(wt_expr + 1)
        interaction_mirnas = log_expr.loc[mirnas]

        fig, ax = plt.subplots(figsize=(1.7, 1.7))

        if params['plot_non_heaps']:
            for val in unfiltered_interaction_df.loc[unfiltered_interaction_df['Gene name'] == wildcards['gene'], 'WT miRNA expression']:
                ax.axvline(np.log2(val + 1), color='#D1e8ff')

        for mirna in mirnas:
            ax.axvline(log_expr.loc[mirna], color='darkblue')
            
        sns.kdeplot(log_expr, ax=ax, color='gray')

        ax.set_xticks([0, max(log_expr)])
        ax.set_xticklabels([0, int(max(wt_expr))])
        ax.set_yticklabels([])
        plt.grid(None)
        sns.despine()
        ax.set_ylabel('distribution')
        ax.set_xlabel('expression')
        # ax.set_title('miRNA expression')
        ax.grid(None)
        plt.tight_layout()
        fig.savefig(output[0])


rule plot_target_log2fc:
    input:
        mrna_data='output/mrna_data_all.csv'
    params:
        vmin=-1.5,
        vmax=1.5
    output:
        'plot/target_data/target_log2fc_{gene}.pdf',
    run:
        # load data
        import matplotlib as mpl
        import seaborn as sns
        import pandas as pd
        import matplotlib.pyplot as plt


        cmap = sns.color_palette("coolwarm", as_cmap=True)
        norm = mpl.colors.Normalize(**params)

        mrna_df = pd.read_csv(input['mrna_data'], header=[0, 1], index_col=[0, 1])
        log2fc = mrna_df.xs('log2FoldChange', axis=1, level=1).droplevel(level=0)[['Dgcr8', 'Drosha', 'Dicer', 'Ago12']]
        log2fc.columns = ['Dgcr8_KO', 'Drosha_KO', 'Dicer_KO', 'Ago2&1_KO']

        gene_data = log2fc.loc[wildcards['gene']]
        color_list = cmap(norm(gene_data))
        fig, ax = plt.subplots(figsize=(2, 2))

        gene_data.plot.bar(color=color_list, ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        sns.despine()
        plt.tight_layout()
        fig.savefig(output[0])

rule plot_target_context_score:
    input:
        interaction_data='output/mirnas/interaction_ranking_all.csv',
        # targetscan_scores='output/targetscan_pairs.csv'
        unfiltered_interaction_data='output/mirnas/unfiltered_interaction_ranking_all.csv'
    output:
        'plot/target_data/ts_scores_{gene}.pdf'
    params:
        mirna_threshold=config['mirna_threshold'],
        plot_non_heaps=True,
    # conda:
    #     "../env/python.yaml"
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        interaction_df = pd.read_csv(input['interaction_data'])
        interactions = interaction_df.loc[interaction_df['Gene name'] == wildcards['gene'], 'context++ score']
        unfiltered_interaction_df = pd.read_csv(input['unfiltered_interaction_data']).query(f'`WT miRNA expression` > {params.mirna_threshold}')
        # ts_scores = pd.read_csv(input['targetscan_scores'])

        fig, ax = plt.subplots(figsize=(1.7, 1.7))


        if params['plot_non_heaps']:
            for val in unfiltered_interaction_df.loc[unfiltered_interaction_df['Gene name'] == wildcards['gene'], 'context++ score']:
                ax.axvline(val, color='#D1e8ff')

        for interaction in interactions:
            ax.axvline(interaction, color='darkblue')

        sns.kdeplot(interaction_df['context++ score'], ax=ax, color='gray')  # TODO maybe use ts_scores

        # ax.set_yscale('log')
        ax.set_xlim(ax.get_xlim()[::-1])
        ax.set_xticks([0, interaction_df['context++ score'].min()]) # ts_scores.min()
        ax.set_yticklabels([])
        plt.grid(None)
        sns.despine()
        ax.set_ylabel('distribution')
        ax.set_xlabel('context++ score')
        # ax.set_title('TargetScan score')
        ax.grid(None)
        plt.tight_layout()
        fig.savefig(output[0])

rule compute_feature_combination_count:
    '''
    How many targets do we get by combining which feature?
    '''
    input:
        raw_interactions='output/mirnas/raw_mirna_interactions_all.csv',
        up_genes='output/up_genes_all.txt',
        mrna_data='output/mrna_data_all.csv'
    output:
        counts='output/feature_combination_count.csv',
        interactions='output/TableS3_Integrative_analysis.xlsx'
    params:
        min_mirna_expression=config['mirna_threshold'],
        # min_mrna_expression=config['mrna_threshold'],
        max_ts_score=config['ts_threshold'],
    conda: '../env/python.yaml'
    script: '../script/compute_feature_combination_count.py'

rule plot_feature_combination_count:
    '''
    How many targets do we get by combining which feature?
    '''
    input:
        'output/feature_combination_count.csv'
    output:
        'plot/feature_combination_count{direction}.svg'
    conda: '../env/python.yaml'
    script: '../script/plot_feature_combination_count.py'

rule plot_interaction_scores:

rule plot_tf_overview:
    input:
        supp_table='output/cluster_kos/SuppTable_mir290_targets.xlsx',
        interaction_ranking='output/mirnas/interaction_ranking_all_istf.csv'
    output:
        'plot/mir_290_regulated_tfs.svg'
    conda: '../env/python.yaml'
    script: '../script/plot_tf_overview.py'

