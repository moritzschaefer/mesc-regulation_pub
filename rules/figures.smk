rule figure1:
    input:
        'plot/ma_all.svg',
        'plot/mrna_pca.svg',
        'plot/feature_combination_count_horizontal.pdf',
        'plot/rip_input_comparison.svg',
        'plot/rip_ago1_ago2_comparison.svg',  # not in manuscript
        'plot/overlap_all_up.svg',
        'plot/overlap_all_down.svg'
        # 'plot/transcript_correlation.svg',
        # 'plot/gene_upregulation.svg',
        # 'plot/gene_downregulation.svg'

rule figure2:
    input:
        expand('plot/target_data/{plot}_{gene}.pdf', plot=['track', 'target_log2fc'], gene=['Prr13', 'Axin2', 'Ctcf', 'Apoe']),
        'plot/interaction_count_histo.svg',
        'plot/interaction_count_gene_expr.svg',
        'plot/interaction_count_mrna_sections.svg',
        'plot/interaction_count_mre_type.svg',
        'plot/heap_for_num_mres.svg',
        'plot/interaction_count_ribo.svg',
        'plot/ribo_seq_validation_all.svg',
        'plot/ms_validation_all.svg',
        'plot/clip_plots/hesc_clip_conservation.svg',
        'plot/clip_plots/hesc_interaction_conservation.svg',
        'plot/mirtarbase_overlap.svg'

rule figure3:
    input:
        'plot/target_cluster_distribution_barplot_all.svg',
        'plot/rip_mir290_loading.svg',
        'plot/cluster_kos/ma_miR-290-295_all.svg',
        'plot/cluster_kos/target_cdf_miR-290-295.svg',
        'plot/cluster_kos/low_up_cdf_supp_miR-290-295.svg',
        'plot/interaction_count_mir290.svg',

rule figure4:
    input:
        'plot/sipool/target_rescue.svg',
        'plot/mir_290_regulated_tfs.svg',
        'plot/hesc_clip_tfap4_track.pdf',
        'plot/target_data/target_log2fc_Tfap4.pdf',  # TODO include into fig 4?
        'plot/mesc_tfap4_track.pdf',
