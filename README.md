# Colorectal-Cancer Early onset 
---
* [Introduction](#Introduction)
* [Resource](#Resource)
* [Pipeline](#Pipeline)

<a name="Introduction"/>

# Introduction
## Identify TF and TF co-oppassoicated with Breast cancer risk
Genome-wide association studies (GWAS) and fine-mapping studies suggest that pathogenic dysregulations of gene expression are mediated through genetically altered binding affinities of specific transcription factors (TFs). Here, we developed a statistical framework to analyze extensive ChIP-seq and GWAS data and identified 22 breast cancer risk-associated TFs. We found that genetic variations affect TF-DNA bindings and the interaction of FOXA1 with co-factors such as ESR1 and E2F1, and the interaction of TFs with chromatin features (i.e., enhancers) played a key role in breast cancer susceptibility. 


## Improved TWAS using the putative regulatory genetic variants
We built gene-expression prediction models using only the putative regulatory genetic variants (flanking 1Mb region) located in the binding sites of these risk-associated TFs at p < 0.05 reported by BCAC (n = 30k SNPs) based on the transcriptome data from GTEx. though we only used 30k SNPs occupied by the 22 selected TFs which shown nominally significant association with breast cancer risk, we were able to predict gene expression at R2>0.01 for 7,538 genes, which is only slightly less than the total number (n=9,109) of predicted genes using all genetic variants in our previous study.  We further focused only on genes that can be predicted by the same set of local genetic variants from either of TCGA and the Molecular Taxonomy of Breast Cancer International Consortium (METABRIC) at R2 > 0.01. we identified 76 genes with predicted expressions that were associated with breast cancer risk at P < 5 × 10-6, at a Bonferroni-corrected significance level as used in our previous study in which we identified 48 genes with regular TWAS. Specifically, we identified 22 genes that are located in regions not yet identified by GWAS (1Mb away). In addition, we uncovered 23 putative breast cancer risk genes in known GWAS loci that had not been previously reported. This work is currently under peer review for publication.

# Resource

### R1. Comprehensive Transcription factor ChIP-seq data in breast cancer cell lines (n = 113) from the ENCODE and the Cistrome database(http://cistrome.org/)

### R2. Summary statistics of imputed SNP data of GWAS from BCAC

### R3.  Gene expression data from TCGA, METABRIC and TCGA

### R4. CRISP-Cas experimental data

### R5. Knock-down experiment data for FOXA1, ESR1, and GATA3


<a name="Pipeline"/>

# Pipeline 
---

## Overview of the Developed Statistical Framework
To investigate how genetic variations of TF-DNA bindings affect breast cancer susceptibility, we developed an analytic framework to analyze ChIP-seq and breast cancer GWAS summary statistics data (Fig. 1A-C). We first identified TF-DNA binding regions by analyzing a total of 113 TF ChIP-seq datasets from multiple breast cancer cell lines collected from ENCODE and the Cistrome database (http://cistrome.org/) (Fig. 1A, B; Table S1; EXPERIMENTAL PROCEDURES). An n × m matrix for n = 11,337,849 genetic variants from the BCAC GWAS data was generated with annotation from m = 113 TF-DNA binding regions. We used Chi-squared values for each genetic variant reported in the BCAC GWAS summary data to measure the association with breast cancer risk. We then used generalized mixed models to estimate the associations between the Chi-squared values (Y) and TF binding status of genetic variants located in binding sites of each TF, given LD blocks of genetic variants to handle the dependence between genetic variants (Fig.1C; Equation 1). To define approximate independent LD blocks similar to other studies (Loh et al., 2015; Pickrell, 2014), we defined LD blocks using non-overlapping segments of 100kb (a similar result with 500kb; see EXPERIMENTAL PROCEDURES).

```
Y_ij=β_0+β_1 〖TF〗_ij+V_i+ε_ij (Equation 1)

```

Specifically, Y_ij is the Chi-square value for the j-th variant in the i-th LD block; β0 is the fixed intercept, and β1 is the fixed slope, which measure the mean difference of the Chi-Square values (∆χ ̅^2) between TF status; 〖TF〗_ij is the j-th TF value (i.e., 1 for a variant located in a TF binding site, 0 otherwise) in the i-th LD block; V_i is the random intercept for the i-th LD block; and ε_ijis the error term. To alleviate the effect of variant correlation within LD blocks on parameter estimation, we divided variants in the whole genome into 100 datasets by sequentially selecting 1 in every 100 variants, according to chromosome positions. We combined results from 100 repeated models using the inverse variance-weighting average method.

## Genetic Variations of TF-DNA Bindings of Breast Cancer Risk-Associated TFs 

## Motif-Dependent Genetic Variations of TF-DNA Bindings of Breast Cancer Risk-Associated TFs

## Genetic Variations of TF-DNA Bindings of FOXA1 and Co-Factors Driving Breast Cancer Susceptibility 

## Genetic Variations of TF Colocalizing with Chromatin Features Associated with Breast Cancer Risk

more information see main.R 
