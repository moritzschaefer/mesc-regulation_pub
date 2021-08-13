rule sipool_pca:
    input:
        'output/sipool/read_counts.tsv'
    output:
        plot='plot/sipool/pca_{subset}.{ext}',
    params:
        sample_colors=config['sample_colors'],
    conda:
        '../env/python.yaml'
    script:
        '../script/cluster_kos/pca.py'

rule sipool_ma_plot:
    input:
        data='output/sipool/TableS7_sipool-Tfap4-Quant-seq.xlsx'
    output: 'plot/sipool/ma_{group}_{subset}.{ext}'
    params:
        baseMean=config['basemean_threshold'],
        padj_threshold=config['padj_threshold']
    conda: '../env/python.yaml'
    script: '../script/cluster_kos/ma_plot.py'

rule sipool_tf_heatmap:
    input:
        quant_seq='output/sipool/read_counts.tsv',
        cluster_ko='output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx'
    output:
        plot='plot/sipool/tf_heatmap.{ext}',
    params:
        genes=list(config['siPOOL_targets'].keys())
    conda:
        '../env/python.yaml'
    script:
        '../script/sipool/tf_heatmap.py'

rule sipool_target_rescue:
    input:
        cluster_ko='output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
        sipool_des='output/sipool/TableS7_sipool-Tfap4-Quant-seq.xlsx',
        interactions='output/mirnas/interaction_ranking_all.csv'
    output:
        plot_png='plot/sipool/target_rescue.png',
        plot_svg='plot/sipool/target_rescue.svg',
        rescue_df='output/sipool/rescues.csv'
    params:
        genes=list(config['siPOOL_targets'].keys()),
        acceptance_factor=2
    conda:
        '../env/python.yaml'
    script:
        '../script/sipool/target_rescue.py'

rule sipool_supp_table:
    input:
        diffexps=lambda wildcards: [str(remote_file(v.format(shrunk=''))) for v in config['siPOOL_targets'].values()],
        expr=remote_file(config['quant_seq_expression_path'])
    output:
        supp_table='output/sipool/TableS7_sipool-Tfap4-Quant-seq.xlsx',
        read_count_matrix='output/sipool/read_counts.tsv'
    params:
        samples=[f'si_{s.lower()}' for s in config['siPOOL_targets'].keys()],
        neg_control='si_neg',
        title='Supp. Table 7: QuantSeq NGS CPM values and differential expression for siPOOL-knockdown of several transcription factors. Referenced by (Supp.) Figure 4'
    conda: '../env/python.yaml'
    script: '../script/supp_table.py'
