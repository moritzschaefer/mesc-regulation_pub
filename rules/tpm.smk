# rule mean_fragment_length:
#     input:
#         bam=config['bam_path']
#     output:
#         txt="stats/{sample}.isize.txt",
#         pdf="stats/{sample}.isize.pdf"
#     log:
#         "log/picard/insert_size/{sample}.log"
#     params:
#         # optional parameters (e.g. relax checks as below)
#         "VALIDATION_STRINGENCY=LENIENT "
#         "METRIC_ACCUMULATION_LEVEL=null "
#         "METRIC_ACCUMULATION_LEVEL=SAMPLE"
#     wrapper:
#         "0.34.0/bio/picard/collectinsertsizemetrics"


rule tpm_normalize_counts:
    input:
        expression='output/rnai_read_counts.tsv',
        annotation="ref/gencode.db",
        # mean_frag_lengths=expand("stats/{sample}.isize.txt", sample=SAMPLES_w2i)
    output:
        "output/all.tpm.tsv"
    conda:
        "../env/python.yaml"
    script:
        "../script/tpm-counts.py"
