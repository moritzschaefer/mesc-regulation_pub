import gffutils

db = gffutils.create_db(snakemake.input[0],
                        dbfn=snakemake.output[0],
                        force=True,
                        keep_order=True,
                        merge_strategy='merge',
                        sort_attribute_values=True)
