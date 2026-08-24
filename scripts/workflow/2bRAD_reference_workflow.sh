#!/usr/bin/env bash

###############################################################################
# Reference-based 2b-RAD population-genomic workflow
# Project: DR_ACER_PopGen_2025
###############################################################################

set -euo pipefail

PROJECT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
cd "$PROJECT"

###############################################################################
# 1. INPUT DATA
###############################################################################

# Raw paired-end 2b-RAD reads:
#   raw_by_sample/N99_R1.fastq.gz
#   raw_by_sample/N99_R2.fastq.gz
#   ...
#   raw_by_sample/N143_R1.fastq.gz
#   raw_by_sample/N143_R2.fastq.gz
#
# Raw reads are deposited separately in NCBI SRA.

###############################################################################
# 2A. MERGE PAIRED-END READS WITH PEAR
###############################################################################

mkdir -p merged_fastq

for r1 in raw_by_sample/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="raw_by_sample/${sample}_R2.fastq.gz"

    pear \
        -f "$r1" \
        -r "$r2" \
        -o "merged_fastq/${sample}"
done

# The assembled reads are used as single-end 2b-RAD fragments downstream.

###############################################################################
# 2B. PREPARE PEAR-ASSEMBLED READS FOR CUTADAPT
###############################################################################

# Decompress PEAR assembled files before Cutadapt if needed.

shopt -s nullglob
assembled_gz=(merged_fastq/*.assembled.fastq.gz)

if (( ${#assembled_gz[@]} > 0 )); then
    for f in "${assembled_gz[@]}"; do
        gzip -cd "$f" > "${f%.gz}"
    done
fi

shopt -u nullglob

###############################################################################
# 3. QUALITY TRIMMING WITH CUTADAPT
###############################################################################

mkdir -p github_shared_trim
mkdir -p logs/cutadapt_shared

for f in merged_fastq/*.assembled.fastq; do
    sample=$(basename "$f" .assembled.fastq)

    python3 -m cutadapt \
        -q 15,15 \
        -m 25 \
        -o "github_shared_trim/${sample}.trim" \
        "$f" \
        > "logs/cutadapt_shared/${sample}_cutadapt.log"
done

###############################################################################
# 4. REFERENCE GENOME PREPARATION
###############################################################################

# Reference:
# Acropora cervicornis jaAcrCerv1.1
# NCBI accession: GCA_964034985.1
#
# Reference FASTA used during the analysis:
GENOME_FASTA="$PROJECT/reference/Acropora_cervicornis.fa"

bowtie2-build "$GENOME_FASTA" "$GENOME_FASTA"
samtools faidx "$GENOME_FASTA"

###############################################################################
# 5. BOWTIE2 ALIGNMENT
###############################################################################

cd "$PROJECT/github_shared_trim"

for f in *.trim; do
    bowtie2 \
        --no-unal \
        --score-min L,16,1 \
        --local \
        -L 16 \
        -x "$GENOME_FASTA" \
        -U "$f" \
        -S "${f}.bt2.sam"
done

###############################################################################
# 6. SORT AND INDEX ALIGNMENTS
###############################################################################

for sam in *.trim.bt2.sam; do
    bam="${sam%.sam}.bam"

    samtools sort \
        -O bam \
        -o "$bam" \
        "$sam"

    samtools index "$bam"
done

###############################################################################
# 7. BUILD BAM LIST
###############################################################################

ls *.trim.bt2.bam |
    sort |
    sed "s|^|$PROJECT/github_shared_trim/|" \
    > "$PROJECT/reference_based/bams_ref"

cd "$PROJECT"

###############################################################################
# 8. ANGSD SNP DISCOVERY AND IBS ANALYSIS
###############################################################################

mkdir -p reference_based/angsd_main_ref

angsd \
    -b reference_based/bams_ref \
    -GL 1 \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -minMapQ 30 \
    -minQ 25 \
    -dosnpstat 1 \
    -doHWE 1 \
    -sb_pval 1e-5 \
    -hetbias_pval 1e-5 \
    -skipTriallelic 1 \
    -minInd 34 \
    -snp_pval 1e-5 \
    -minMaf 0.05 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -doCounts 1 \
    -makeMatrix 1 \
    -doIBS 1 \
    -doCov 1 \
    -doGeno 8 \
    -doBCF 1 \
    -doPost 1 \
    -doGlf 2 \
    -P 1 \
    -out reference_based/angsd_main_ref/ref_main

# Save retained SNP coordinates and alleles for downstream analyses.

zcat reference_based/angsd_main_ref/ref_main.mafs.gz |
    awk 'NR > 1 {print $1, $2, $3, $4}' OFS='\t' \
    > reference_based/retained_2515_sites_major_minor.txt

n_sites=$(wc -l < reference_based/retained_2515_sites_major_minor.txt)

if [ "$n_sites" -ne 2515 ]; then
    echo "Expected 2515 retained SNPs but found $n_sites" >&2
    exit 1
fi

###############################################################################
# 9. IBS PCoA
###############################################################################

Rscript publication_repository/scripts/analysis/angsd_ibs_pca_fixed.R \
    reference_based/angsd_main_ref/ref_main.ibsMat \
    reference_based/bams_ref \
    "" \
    reference_based/ref_main

# The first two PCoA axes explain:
#   PCoA1 = 46.85607%
#   PCoA2 = 29.77815%

###############################################################################
# 10. GENET ASSIGNMENT
###############################################################################

Rscript publication_repository/scripts/analysis/assign_genets.R \
    reference_based/angsd_main_ref/ref_main.ibsMat \
    reference_based/bams_ref \
    reference_based/genotype_assignments.tsv

###############################################################################
# 11. FULL 45-COLONY ngsRelate ANALYSIS
###############################################################################

cd "$PROJECT/reference_based/angsd_main_ref"

# IDs must follow the same ordering as the BAM/BCF order.
cut -f1 "$PROJECT/reference_based/bams_ref" |
    sed 's#.*/##' |
    sed 's/\.bam$//' \
    > ref_ids.txt

ngsRelate \
    -h ref_main.bcf \
    -n 45 \
    -z ref_ids.txt \
    -O ref_kinship

cd "$PROJECT"

###############################################################################
# 12. REDUCED ONE-REPRESENTATIVE-PER-GENET ANALYSIS
###############################################################################

# Representative ramets:
#
#   N120 = AC3
#   N128 = AC2
#   N130 = AC1
#
# Reduced analysis using one representative per genet.

cat > reference_based/reduced_3genets.bams.txt <<BAMS
$PROJECT/github_shared_trim/N120.trim.bt2.bam
$PROJECT/github_shared_trim/N128.trim.bt2.bam
$PROJECT/github_shared_trim/N130.trim.bt2.bam
BAMS

mkdir -p reference_based/angsd_reduced_3genets

angsd \
    -b reference_based/reduced_3genets.bams.txt \
    -GL 1 \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -minMapQ 30 \
    -minQ 25 \
    -sites reference_based/retained_2515_sites_major_minor.txt \
    -doMajorMinor 3 \
    -doMaf 1 \
    -doCounts 1 \
    -doPost 1 \
    -doGeno 8 \
    -doBCF 1 \
    -doGlf 2 \
    -minInd 1 \
    -P 1 \
    -out reference_based/angsd_reduced_3genets/reduced

# Sample IDs in BAM/BCF order.
cut -f1 reference_based/reduced_3genets.bams.txt |
    sed 's#.*/##' |
    sed 's/\.bam$//' \
    > reference_based/angsd_reduced_3genets/reduced_ids_bcf_order.txt

ngsRelate \
    -h reference_based/angsd_reduced_3genets/reduced.bcf \
    -n 3 \
    -z reference_based/angsd_reduced_3genets/reduced_ids_bcf_order.txt \
    -O reference_based/ngsrelate_reduced_3genets.txt

###############################################################################
# END
###############################################################################
