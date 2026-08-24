# Genomic characterization of *Acropora cervicornis* nurseries in the Dominican Republic

This repository contains analysis scripts, processed genomic data, metadata, and supplementary outputs associated with the study:

**Genomic characterization reveals high clonal redundancy in two *Acropora cervicornis* nurseries in the Dominican Republic**

Shamwari Anseeuw Carrasco, Kasey Walsh, Rebecca Garcia-Camps, Ainhoa L. Zubillaga, Aldo Croquer, and Debashish Bhattacharya

## Overview

We used 2b-RAD sequencing to characterize multilocus genotypic variation among 45 *Acropora cervicornis* colonies maintained in two in-situ coral nurseries in Punta Cana, Dominican Republic.

Reference-based SNP discovery and filtering produced a final dataset of 2,515 SNPs. The analyses included principal coordinates analysis (PCoA), average-linkage hierarchical clustering, pairwise identity-by-state (IBS) distances, genotypic diversity, and relatedness estimation.

## Data availability

Raw 2b-RAD sequencing reads are available through the NCBI Sequence Read Archive under BioProject **PRJNA1514087**.

Raw FASTQ files are not included in this repository. Processed data and analysis scripts required to reproduce the analyses presented in the manuscript are provided here.

## Repository structure

```text
.
├── metadata/          Sample metadata and genet assignments
├── processed_data/    Processed genomic data and analysis outputs
├── scripts/
│   ├── analysis/      Core analysis scripts
│   ├── publication/   Figure and supplementary-table scripts
│   ├── validation/    Analysis validation scripts
│   └── workflow/      Reference-based 2b-RAD workflow
├── figures/           Manuscript figures
└── supplementary/     Supplementary tables and figures
```

## Analysis workflow

The reference-based workflow is documented in:

```text
scripts/workflow/2bRAD_reference_workflow.sh
```

The workflow includes read merging and trimming, alignment to the *Acropora cervicornis* jaAcrCerv1.1 reference genome (NCBI accession **GCA_964034985.1**), SNP discovery and IBS analysis with ANGSD, PCoA, hierarchical clustering and genet assignment, and relatedness analysis with ngsRelate.

Core downstream analyses are contained in `scripts/analysis/`.

## Reproducing figures and supplementary tables

Scripts used to generate publication outputs are contained in:

```text
scripts/publication/
```

From the repository root:

```bash
Rscript scripts/publication/Figure2.R
Rscript scripts/publication/Figure3.R
Rscript scripts/publication/FigureS1.R

Rscript scripts/publication/TableS1.R
Rscript scripts/publication/TableS2.R
Rscript scripts/publication/TableS3.R
Rscript scripts/publication/TableS4.R
```

These scripts use the metadata and processed data included in this repository and do not require the raw FASTQ files.

## Citation

Citation information will be added following publication.
