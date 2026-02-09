#!/usr/bin/env bash
set -euo pipefail


THREADS=48


############################################
# 01. Genes in Highly Divergent Regions
############################################
minimap2 -c -x asm5 -o hap1_hap2.paf -t ${THREADS} --eqx hap1.fa hap2.fa
syri -F P -c hap1_hap2.paf -r hap1.fa -q hap2.fa --nc 12 --prefix hap1_hap2


grep HDR hap1_hap2syri.out | cut -f 1,2,3,9 > hdr.bed
bedtools intersect -a hap1_gene.bed -b hdr.bed > candidate_gene.bed
cut -f 4 candidate_gene.bed > candidate_gene.id


############################################
# 02. Genes with Specific Expression
############################################
awk 'NR==1{print "Geneid\tLeaf\tStigma\tPollen"; next}
{leaf=($2+$3)/2; stigma=($4+$5)/2; pollen=($6+$7)/2;
print $1"\t"leaf"\t"stigma"\t"pollen}' all_tpm > all_tpm_avg.tsv


grep -F -f candidate_gene.id all_tpm_avg.tsv > candidate_gene.TPM


awk '$2 < 1 && $3 >= 1 && $4 < 1' candidate_gene.TPM > candidate_gene_pistil.TPM
awk '$2 < 1 && $3 < 1 && $4 >= 1' candidate_gene.TPM > candidate_gene_pollen.TPM


############################################
# 03. Genes with Linkage
############################################
cut -f 1 candidate_gene_pollen.TPM > HDR_pollen.id
cut -f 1 candidate_gene_pistil.TPM > HDR_style.id


grep -f HDR_pollen.id hap1_gene.bed > HDR_pollen.bed
grep -f HDR_style.id hap1_gene.bed > HDR_style.bed


awk '{print $0"\tPollen"}' HDR_pollen.bed > HDR_pollen_anno.bed
awk '{print $0"\tStyle"}' HDR_style.bed > HDR_style_anno.bed


cat HDR_pollen_anno.bed HDR_style_anno.bed | sort -k1,1 -k2,2n > HDR_anno.bed


bedtools merge -i hap1_gene.bed > gene.merged.bed
cut -f1,2 hap1.fa.fai > genome.size
bedtools complement -i gene.merged.bed -g genome.size > intergenic.bed


awk '{sum += $3 - $2; n++} END {print "Average intergenic length:", sum/n}' intergenic.bed


# The average intergenic length should be added to find_link.py
python3 find_link.py HDR_anno.bed > clusters.txt
cut -f 4 clusters.txt | sort | uniq > link.geneid


############################################
# 04. Genes with Subcellular Location
############################################
biolib run DTU/DeepTMHMM --fasta hap1_pep.fa
# Manual filtering of localization results is required


############################################
# 05. Allelic Gene Diversity
############################################
blastp -query hap1_can_pep.fa -db hap2_pep.fa \
-outfmt 6 -out hap1_can.blastp -num_threads 4
# Manual pairing and extraction of allelic gene pairs is required
