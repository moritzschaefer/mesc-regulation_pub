# README

## Install & import

To facilitate running dedicated parts of the pipeline to generate specific plots, corresponding snakemake commands and script names are provided on a per-figure basis (including supplementary figures). All you need to run the pipeline is conda (4.9.1), snakemake (5.26.1) and my support library. Running the following commands should get you to a running environment:

You can initialize the environment represented by the conda.yaml file.

``` shell
conda env create -f env/python.yaml --name mesc
```

Next, initialize pyensembl and clone the repository (if you haven't done so yet)
``` shell
pyensembl install --release 98 --species mus_musculus
cd ~/  # or wherever you want to store the source code 
git clone https://github.com/moritzschaefer/mesc-regulation_pub
cd mesc-regulation
```


By running `snakemake -j 4 --snakefile import_supptables.snakefile` you can import Supplementary Tables as provided by the applications. This way, there is no need to rerun any NGS pipelines on raw (fastq) data files. Make sure to place all supplementary tables as well as the quantseq-read counts file from GEO in the `supp/` folder.

Note that the generation of all output files might take a substantial amount of time and RAM!

## Run 

Running the pipeline is as simple as the following (here for generating all plots from Figure 1)

``` shell
snakemake -j 4 figure1
```

## Additional information

- To generate the genome track plots from figure 2, please manually fill the `file =` paths in misc/tracks_manual.ini
- The TPM calculation *ignores* mean fragment length, because the resulting TPM read counts wouldn't correspond to log2FCs from DESeq2 anymore. This is somewhat inaccurate but it is the most practical solution.
- snakePipes' mRNA-seq pipeline was used to generate read count tables and DEG analyses
