<img align="center" src="docs/source/images/hich_logo_white.png" alt="Hich" width="100" height="100" style="margin-right: 10px;">
Hich is a powerful workflow for 3D Genome (Hi-C) data. It replaces previous tools such as Juicer, HiCPro, nf-core hic, and distiller.

# [Documentation and vignettes](https://hich.readthedocs.io/en/latest/index.html)

Hich installs with conda/mamba. Sensible defaults are chosen and can be easily modified with on a granular, per-sample level for total control. Setup for an experiment is based on one sample file that looks like this:

| | | | | | | |
|-|-|-|-|-|-|-|
|condition|biorep|techrep|fastq1|fastq2|assembly|enzymes|
|H|1|1|Human_1_1_R1.fq.gz|Human_1_1_R2.fq.gz|hg38|Arima|
|H|1|2|Human_1_2_R1.fq.gz|Human_1_2_R2.fq.gz|hg38|Arima|
|H|2|1|Human_2_1_R1.fq.gz|Human_2_1_R2.fq.gz|hg38|Arima|
|H|2|2|Human_2_2_R1.fq.gz|Human_2_2_R2.fq.gz|hg38|Arima|
|M|1|1|Mouse_1_1_R1.fq.gz|Mouse_1_1_R2.fq.gz|mm10|Arima|
|M|1|2|Mouse_1_2_R1.fq.gz|Mouse_1_2_R2.fq.gz|mm10|Arima|
|M|2|1|Mouse_2_1_R1.fq.gz|Mouse_2_1_R2.fq.gz|mm10|Arima|
|M|2|2|Mouse_2_2_R1.fq.gz|Mouse_2_2_R2.fq.gz|mm10|Arima|



# Features

✔ **Easy.** Ultra portable and installs in seconds. Experiments can be set up in a couple minutes. Hich allows stub runs and can also auto-downsample your input data for a quick "humid" test before a full run. The run can be paused and resumed after any step, and successfull intermediate outputs are cached, so that common issues in an HPC environment such as running out of time or memory can be rerun from where they left off. Hich runs in parallel to take full advantage of the power of modern HPC environments, but it can also be run on a personal laptop.

✔ **Comprehensive.** Produce compartment scores, differential loop analysis, TADs and contact matrices from raw paired-end sequencing reads (.fastq) or intermediate file formats (.sam, .bam, .pairs). Unique downsampling feature matches coverage for differential feature calling to ensure results are not biased by sequencing depth. Features can be called with multiple parameters for sensitivity analysis, and sensible defaults are chosen.

✔ **Replicates and controls.** Manage multi-condition experiments by merging multiple technical and biological replicates

✔ **Browseable.** Outputs compatible with both popular Hi-C browsers, [Juicebox](https://www.aidenlab.org/juicebox/) and [higlass](https://higlass.io/)

✔ **Filters.** Eliminates PCR duplicates, undigested chromatin, and other technical artifacts.

✔ **Compatible with advanced assays and commercial kits.** Works with conventional Hi-C, capture-based methods, single- and multi-restriction enzyme digests, and non-restriction enzyme digests such as MNase and DNase. Specialized settings for commercial kits from [Arima](https://arimagenomics.com/products/genome-wide-hic/), [Dovetail](https://cantatabio.com/dovetail-genomics/products/), and [Phase](https://phasegenomics.com/products/proximo/). Can also align Hi-C reads that have been subjected to deamination via bisulfite or enzymatic conversion to study the methylome.

✔ **Advanced QC.** Novel Hicrep-based clustermapping for QC to compare conditions and confirm similarity of replicates at chromosome scale 
