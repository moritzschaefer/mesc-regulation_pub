# see also /home/moritz/SparkleShare/github.com/wiki/roam/20210209073259-mir_290ko_sequencing_analysis.org

rule cluster_kos_pca:
    input:
        'output/cluster_kos/read_counts.tsv'
    output:
        plot='plot/cluster_kos/pca_{subset}.{ext}',
    params:
        sample_colors=config['sample_colors'],
    conda:
        '../env/python.yaml'
    script:
        '../script/cluster_kos/pca.py'

rule cluster_ma_plot:
    input:
        data='output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx'
    output: 'plot/cluster_kos/ma_{group}_{subset}.{ext}'
    params:
        baseMean=config['basemean_threshold'],
        padj_threshold=config['padj_threshold']
    conda: '../env/python.yaml'
    script: '../script/cluster_kos/ma_plot.py'

rule cluster_target_cdf:
    input:
        data='output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
        mirna_targets='output/mirnas/interaction_ranking_all.csv',
        low_mirna_targets='output/mirnas/low_interaction_ranking_all.csv',
        up_genes='output/up_genes_all.txt',
        low_up_genes='output/low_up_genes_all.txt',
        ago2_heap='output/ago2_heap_mesc_mres.csv',
        unfiltered_interaction_data='output/mirnas/unfiltered_interaction_ranking_all.csv',
        mrna_data='output/mrna_data_all.csv'
    output:
        plot='plot/cluster_kos/target_cdf_{cluster}.svg',
        supp_plot='plot/cluster_kos/target_cdf_supp_{cluster}.svg',
        low_up_supp_plot='plot/cluster_kos/low_up_cdf_supp_{cluster}.svg',
        up_percentage='output/cluster_percentage_{cluster}.txt'
    params:
        sample_colors=config['sample_colors'],
        mirna_threshold=config['mirna_threshold'],
        plot_ts_predictions=True,
        plot_heap_predictions=True,
        plot_upregulation_predictions=True
    conda: '../env/python.yaml'
    script: '../script/cluster_kos/target_cdf.py'

clusters = list(config['cluster_mutants'].keys())
rule cluster_ko_supp_table:
    input:
        diffexps=lambda wildcards: [str(remote_file(config['cluster_mutants'][cluster]['file'].format(shrunk='_LFCshrunk'))) for cluster in clusters],  # use LFCshrunk in the miRNA KO quant-seq
        expr=remote_file(config['quant_seq_expression_path'])
    output:
        supp_table='output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
        read_count_matrix='output/cluster_kos/read_counts.tsv'
    params:
        samples=clusters,
        neg_control='WT',
        title='Supp. Table 6: QuantSeq NGS CPM values and differential expression for miRNA cluster KO mutant. Referenced by (Supp.) Figures 3 and 4'
    conda: '../env/python.yaml'
    script: '../script/supp_table.py'

rule cluster_mir290_targets_supp_table:
    input:
        mir_cluster_des='output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
        interactions='output/TableS3_Integrative_analysis.xlsx',
        tf_annotation='output/tf_annotation.csv',
    output:
        supp_table='output/cluster_kos/SuppTable_mir290_targets.xlsx',
    conda: '../env/python.yaml'
    script: '../script/cluster_kos/mir290_targets_supp_table.py'
