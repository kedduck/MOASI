# Multi-Omics Assisted S-locus Identification

This repository provides a reproducible pipeline for identifying candidate **S-locus genes** using multi-omics evidence, including genome structural variation, tissue-specific expression, genetic linkage, subcellular localization, and allelic diversity.

The workflow is designed for **haplotype-resolved assemblies** and is suitable for highly divergent genomic regions (HDRs).

---

## Overview of the Pipeline

1. **Genes in Highly Divergent Regions (HDRs)**  
   Identify genes located in structurally divergent regions between two haplotypes.

2. **Genes with Specific Expression**  
   Filter HDR genes based on tissue-specific expression patterns.

3. **Genes with Genetic Linkage**  
   Cluster candidate genes based on physical proximity and linkage signals.

4. **Genes with Subcellular Location**  
   Predict transmembrane domains to infer protein localization.

5. **Allelic Gene Diversity**  
   Assess allelic divergence between haplotypes using protein sequence similarity.

---

## Requirements

### Software

- minimap2
- syri
- bedtools
- blast+
- awk, grep, cut, sort (GNU coreutils)
- Python ≥ 3.7
- DeepTMHMM (via BioLib)

### Input Files

- `hap1.fa`, `hap2.fa` – haplotype genome assemblies
- `hap1_gene.bed` – gene annotation (haplotype 1)
- `hap1.fa.fai` – genome index (samtools faidx)
- `hap1_pep.fa`, `hap2_pep.fa` – protein sequences
- `all_tpm` – gene expression TPM table

---

## Notes on Manual Steps

- **DeepTMHMM results**: manually inspect transmembrane predictions to select membrane-localized candidates.
- **Allelic diversity analysis**: manually pair haplotype-specific alleles based on BLASTP results and sequence similarity.
- **find_link.py**: the average intergenic length calculated above must be added as a parameter inside the script.
