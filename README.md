# Analysis Scripts for the Microbiota--Host Transcriptomic Study

## Overview

This repository contains the core scripts used for RNA-seq differential
expression analysis, WGCNA, ATAC-seq processing, and alternative
splicing analysis in this study.

The scripts are provided to document the main computational workflows
and facilitate reproducibility. Processed results are reported in the
main text and Supplementary Tables; this repository focuses on the
analysis code rather than duplicating those result files.

## Repository contents

  -----------------------------------------------------------------------
  File                                Description
  ----------------------------------- -----------------------------------
  `01.DESeq2.R`                       Differential expression analysis of
                                      bulk RNA-seq data using DESeq2.

  `02.WGCNA.R`                        Weighted gene co-expression network
                                      analysis, including network
                                      construction, module detection, and
                                      module--trait association analysis.

  `03.ATAC.sh`                        ATAC-seq read-processing workflow.

  `04.rmats.sh`                       rMATS workflow for alternative
                                      splicing analysis.

  `05.rmats_sig_events.R`             Processing and filtering of
                                      significant rMATS alternative
                                      splicing events.
  -----------------------------------------------------------------------

## Tissue-specific analyses

The same analysis scripts were applied independently to the relevant
tissues. In particular, the RNA-seq differential expression and WGCNA
workflows were run separately for each tissue using the corresponding
tissue-specific expression/count matrix and sample/trait information.

The scripts were not run by pooling samples from different tissues into
a single analysis. For each tissue, the input file paths, sample
information, and output directory should be set to the corresponding
tissue-specific files. The analysis code and analytical parameters are
otherwise kept consistent across tissues, unless explicitly stated in
the manuscript.

To reproduce a tissue-specific analysis:

1.  Prepare the input matrix and sample/trait metadata for the tissue of
    interest.
2.  Check that sample identifiers and their order match between the
    expression/count matrix and metadata.
3.  Set the input and output paths in the relevant script.
4.  Run the same script for that tissue.
5.  Repeat the procedure independently for each additional tissue.

The ATAC-seq and rMATS scripts should be run using the corresponding
samples and comparison specified in the manuscript. Do not assume that
these workflows were performed for every tissue unless indicated by the
study design.

## General requirements

The scripts require the relevant software and R packages to be installed
and available in the execution environment. Depending on the workflow,
these include:

-   R and the R packages required by the individual script (for example,
    DESeq2 and WGCNA).
-   The command-line tools required by the ATAC-seq and rMATS workflows.
-   Reference genome files, gene annotation files, and input data
    corresponding to the versions used in the study.

Please check the package and software requirements at the beginning of
each script and install the required dependencies before execution.

## Input data and file paths

Input file names and directory paths may need to be adapted to the local
computing environment. Before running a script, replace any example or
project-specific paths with the paths to the corresponding input files.

For analyses involving multiple tissues, use the matching
tissue-specific input files and metadata for each run. Ensure that
sample IDs are consistent across all input files.

## Outputs

The scripts generate intermediate files, statistical result tables,
and/or plots according to the commands in each script. The principal
processed results used in the manuscript are presented in the article
and Supplementary Tables.

## Data availability

Raw sequencing data are available through the CNGB Sequence Archive
(CNSA) of the China National GeneBank Database (CNGBdb) under accession
numbers **CNP0009911** and **CNP0005541**. Processed datasets and
analytical results are provided in the article and Supplementary
Materials.

## Notes

-   These scripts are shared for research reproducibility.
-   Results may differ if software versions, reference annotations,
    genome builds, input data, or analysis parameters are changed.
-   Please consult the manuscript and Supplementary Methods for the
    experimental design, sample information, and analysis criteria.
