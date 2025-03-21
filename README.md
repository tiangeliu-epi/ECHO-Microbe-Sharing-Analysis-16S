# ECHO-Microbe-Sharing-Analysis-16S

This repository contains the R scripts and RMarkdown files used to analyze mother-child microbe sharing for the article: Liu, Tiange, et al. "Maternal Vaginal and Fecal Microbiota in Later Pregnancy Contribute to Child Fecal Microbiota Development in the ECHO Cohort." iScience (2025). https://www.cell.com/iscience/fulltext/S2589-0042(25)00472-9 

---
**Repository Structure:**
---
The repository is organized into folders representing key analysis stages:

1. Function: Shared functions used across scripts.
2. Data Prep & Table 1: Scripts for preparing data and summarizing participant characteristics.
3. Analysis: Scripts for analyzing data for each study Aim and in each cohort.
---
**Description of Scripts:**
---
1. 0_Function.R

Contains reusable functions for data processing, statistical analyses, and visualization. Need to source this script before running any other scripts.

2. 1_Read_covariates_sql.R

Processes and cleans metadata from SQL databases. The resulting metadata is used for all the analyses.

3. 2.1_Table_Sequencing.Rmd

Summarizes information on microbiota data collection and sequencing techniques in each cohort, as presented in Supplement Table 2.

4. 2.2_Table.Rmd

Generates Table 1 of the manuscript, summarizing participant characteristics.

5. 3.1_MARCH_FEAST.Rmd and 3.2_MARCH_FEAST_summary.Rmd

Conducts microbial source tracking using the FEAST algorithm in the MARCH cohort for Aim 1 and Aim 3. Results are presented in Figure 2 D, Supplement Figure 4 D, and Supplement Figure 6 B1-B3.

6. 3.3_MARCH_overlap.Rmd

Summarizes the number of amplicon sequence variants (ASVs) that overlap between maternal and infant samples in the MARCH cohort for Aim 1. Results are presented in Figure 2 A0-A3, and Supplement Figure 4 A0-A3.

7. 3.4_MARCH_FEAST_randomize.Rmd

Estimates the contribution of biological vs. random mother to infant fecal microbiota in the MARCH cohort, as a sensitivity analysis for Aim 1. Results are presented in Supplement Figure 5.

8. 4_MARCH_ASV.Rmd

Conducts the ASV sharing analysis in the MARCH cohort for Aim 1 and Aim 3. Results are presented in Figure 2 B0-B3 & C, Supplement Figure 4 BO-B3 & C, and Supplement Figure 6 A1-A3. This script also analyzes the sharing of Bifidobacterium and Bacteroides in the MARCH cohort, results of which are presented in Figure 4, Supplement Figure 7, and Supplement Figure 8.

9. 5.1_RESONANCE_FEAST.Rmd, 5.2_ROCHESTER_FEAST.Rmd, 5.3_VDAART_FEAST.Rmd, 5.4_WISC_FEAST.Rmd

Conduct microbial source tracking using the FEAST algorithm in the RESONANCE, Rochester, VDAART, and WISC cohort for Aim 2.

10. 6_Aim3_FEASTsum_ASV.Rmd

Summarizes and visualizes FEAST results obtained in 5.1-5.4 scripts, and conducts ASV sharing analysis in the RESONANCE, Rochester, VDAART, and WISC cohort for Aim 2 and Aim 3. Results are presented in Figure 3. This script also contains analysis on the sharing of Bifidobacterium and Bacteroides in these 4 cohorts, results of which are presented in Figures 5 & 6.

11. 7_Similarity.Rmd

Compares microbiota dissimilarity between child-biological mother and child-random mother dyads across 5 cohorts. Results are presented in Supplement Figure 3.

12. 8_Maternal_abundant_genera.Rmd

To plot the dominant genera of maternal fecal and vaginal microbiota in each cohort. Results are presented in Supplement Figure 2.

---
**Data Availability**
---
Select de-identified data from the ECHO Program are available through NICHD’s Data and Specimen Hub (DASH) (https://dash.nichd.nih.gov/). Information on study data not available on DASH, such as some Indigenous datasets, can be found on the ECHO study DASH webpage (https://dash.nichd.nih.gov/explore/study?q=echo&filters=%5b%5d&page=1&sortBy=relevance&asc=true&size=50). The 16S rRNA amplicon gene sequences are publicly available via accession numbers listed in the Key Resources Table in the published article. 
