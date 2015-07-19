# Systematic Rational Identification of Sex-Linked Molecular Alterations and Therapies in Cancer
## Sadhika Malladi<sup>1</sup>, Jonathan Q. Ma<sup>1</sup>, Dr. Andrew Beck

### Introduction
We have created this repository to house the code relevant to the project "Systematic Rational Identification of Sex-Linked Molecular Alterations and Therapies in Cancer." The repository contains R scripts used to process and anaylyze data and results.

Our objectives were to identify genomic differences between male and female cancer patients and extrapolate these findings to predict which drugs would be efficient in males, females, or both. Currently, sex is recognized as a factor in the incidence and progression of many diseases, but few studies have investigated the differences to aid in the development of personalized treatment. Moreover, many of the drugs in use today were developed during a time when the Food and Drug Administration (FDA) excluded women from phase I and II human clinical trials, suggesting that females may be receiving sub-optimal treatment due to their underrepresentation during testing phases. Even as of 2010, only 37% of trial participants in 9 major medical journals were female, suggesting an enduring trend.

In our work, we identify cancer-specific and pan-cancer genomic differences between male and female patients. Then, we use Ingenuity Pathway Analysis (IPA) to discover pathways variably enriched in male and female tumors, lending insight into the biological processes that differ between both sexes. Finally, we use our genomic findings to identify drugs that work well in males, females, or both sexes.

### Navigation
We have divided our scripts into 6 different folders, each containing a different analysis. If you want to replicate our results, however, you will need the data.

#### Expression Analysis
This folder contains the scripts for foundational analysis of gene expression differences in male and female tumors. The script `Complete_Expression_Analysis.R` contains the p-value computation and adjustment, along with fold change and meta-analysis functionality. The results from this script provide the basis for the rest of the analyses. We also adjust the p-values within tissues separately (instead of across all tissues) in `Adjust_PValues.R`. 

#### Copy Number Analysis
This folder contains the scripts for analyzing male-female differences in copy number within and across cancers. First, `Download_Data.R` downloads data from cBioPortal (an outlet of TCGA). Then, the p-values and fold changes are calculated in the same manner as in `Expression Analysis` in `Complete_Analysis.R`. To analyze the significance of the differences in copy number, `Identify_Significant.R` isolates genes with a p-value less than 0.05. Moreover, the script cross-references the findings from `Expression Analysis` with the copy number findings to see if any genes demonstrate significant differences in RNA expression and DNA copy number. We complement the significance measurements with tissue-specific fold change computations in `TissueSpecific_FoldChanges.R`.

#### Methylation Analysis
This folder contains the scripts for analyzing sex differences in DNA methylation both within and across cancers. 'Download_Cbio_Methyl_Data.R' downloads relevant methylation data from cBioPortal. 'Main_Analysis.R' analyzes this data using the wilcoxon rank-sum test, yielding Benjamini-Hochberg-values which represent the significance of sex disparities in methylation. 'Fisher_analysis.R' aggregates cancer-specific p-values into pan-cancer p-values. 'FoldChanges.R' computes fold changes of methylation with respect to sex. Finally, human_genes.txt can be modified to only perform the analysis on a subset of available HUGO genes.

#### Resistance and Sensitivity Signatures
This folder contains the scripts required to define cancer-specific and pan-cancer resistance and sensitivity signatures for male and female patients separately. `PreProcessing.R` and `TumorvNormal.R` comptue p-values and fold changes summarizing the differences between neoplastic and nonneoplastic tissue in males and females separately. Thus, the scripts yield resistance and sensitivity signatures for each sex in each tissue. These signatures are used as input to IPA and LincsCloud. An additional script, `IngenuityInput.R`, formats the signatures for use with IPA. The LincsCloud inputs were simply copied and pasted from the generated text files into the web service.

#### Interpret Results
This folder contains the scripts required to process and analyze the results from IPA and LincsCloud. The scripts generate useful plots and tables that assist in finding patterns in the results.

#### Permutation Analysis
Performs a permutation-based analysis to validate LincsCloud results. 'Create_Permuted_Gene_Signature.R' essentially re-executes our Gene Expression Analysis with sex permuted at the sample level and permuted p-values computed using the Wilcoxon Rank-Sum Test through Significance Analysis of Microarrays, in order to yield permuted signatures under the null hypothesis that sex does not affect drug efficacy. Next, 'permutedqueryscript_SAM.csh' and 'extraneousdelete.csh' execute LincsCloud analysis on these permuted signatures. These raw results are then interpreted by 'ResultAnalyzer.R' to yield empirical Benjamini-Hochberg values, quantifying the significance of observed sex disparities in connectivity score.

NOTE: 'permutedqueryscript_SAM.csh' and 'extraneousdelete.csh' must be run on the LINCS Compute Connectivity on the Cloud (C3) interface, accessible through SSH. To use LINCS C3, you must request access by emailing lincscloud@broadinstitute.org. For more information, visit c3.lincscloud.org.