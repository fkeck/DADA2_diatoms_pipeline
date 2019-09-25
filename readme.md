# DADA2 pipeline for rbcL diatoms

## Description

This is a custom script to process Illumina MiSeq HTS data directly modified from the workflows published
on the official [DADA2 website](https://benjjneb.github.io/dada2/index.html). The script includes some modifications to fulfill
specific needs and has proven to work well for diatom metabarcoding with rbcL. Taxonomic assignment step is done using Diat.barcode (Rimet et al. 2018), an expertly curated barcode library for diatoms.

## Dependencies

The pipeline require some packages to be installed.

- `ShortRead` and `dada2` from Bioconductor:

    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("ShortRead")
    BiocManager::install("dada2")

- `ggplot2` from CRAN:

    install.packages("ggplot2")
    
- `diatbarcode` from GitHub:

    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
    devtools::install_github("fkeck/diatbarcode")
    
You also will need a recent version of `cutadapt` if you have to remove your primers (first step of the pipeline). Check the [cutadapt website](https://cutadapt.readthedocs.io/en/stable/installation.html) for installation instructions. Installing cutadapt on Windows seems to be possible but tricky.

## References

> Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods, 13(7), 581.
    
> Rimet, Frederic; Gusev, Evgenuy; Kahlert, Maria; Kelly, Martyn; Kulikovskiy, Maxim; Maltsev, Yevhen; Mann, David; Pfannkuchen, Martin; Trobajo, Rosa; Vasselon, Valentin; Zimmermann, Jonas; Bouchez, Agnès, 2018, "Diat.barcode, an open-access barcode library for diatoms", https://doi.org/10.15454/TOMBYZ, Portail Data Inra, V4
    
> Keck, François; Rimet, Frédéric; Vasselon, Valentin; Bouchez, Agnès, 2019, "A ready-to-use database for DADA2: Diat.barcode_rbcL_312bp_DADA2", https://doi.org/10.15454/HNI1EK, Portail Data Inra, V1 