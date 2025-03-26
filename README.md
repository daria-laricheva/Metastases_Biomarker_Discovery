# Peritoneal Metastases Biomarker Discovery

## Goals
The primary objective of this project is to identify gene biomarkers that distinguish peritoneal metastases from primary colorectal tumors. This identification is crucial as peritoneal metastasis significantly alters the prognosis and treatment options for colorectal cancer patients. Early detection and a deeper understanding of the molecular characteristics of peritoneal metastases could greatly enhance patient management.

## Necessary tools
### Tools utilized to successfully run `pipeline.py`:
1. `FastQC (v0.12.1)`
2. `MultiQC (v1.22.2)`
3. `HISAT2 (v2.2.1)`
4. `samtools (v1.19.2)`
5. `featureCounts (v2.0.6)`, a part of the `Subread` package

### Key R packages used 
(For the complete list of packages and versions used in the R environment please see `R_session_info.txt`):
1. `DESeq2 (v1.38.3)`
2. `ComplexHeatmap (v2.14.0)`
3. `ggplot2 (3.5.1)`

## Files description
1. `pipeline.py` -- python script for RNA-Seq data preprocessing/processing. Includes steps such as:
     * quality control
     * building the index for a reference genome
     * performing the alignment of reads
     * converting a SAM file to a sorted BAM file and indexing it
     * subsetting reference genome file to chr17 (because of memory constraints)
     * gene expression quantification
2. `Metastases_biomarker_discovery.Rmd` -- R Markdown file containing the full differential gene expression analysis
3. `Metastases_biomarker_discovery.pdf` -- PDF rendering of the R Markdown file, including visualizations and results
4. `Metastases_biomarker_discovery-Project_report.pdf` -- comprehensive project report detailing the rationale, step-by-step methodology, tool versions used, and interpretation of results
5. `filtered_significant_genes.csv` -- CSV file containing significant genes identified through differential gene expression analysis
6. `R_session_info.txt` -- TXT file containing the complete list of packages and versions used in the R environment

## Input
### The `pipeline.py` expects:
1. Paired-end RNA-seq FASTQ files
2. Reference genome file
3. Gene annotation file

### The `Metastases_biomarker_discovery.Rmd` expects:
1. Combined counts file

## Output
### The `pipeline.py` outputs:
1. FastQC quality control reports for each input FASTQ file
2. Summary MultiQC report aggregating all FastQC outputs
3. Genome index files folder
4. SAM files generated from aligning RNA-seq reads
5. Sorted and indexed BAM files derived from the aligned SAM files
6. Count matrix for downstream differential expression analysis
