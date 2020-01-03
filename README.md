# Supplementary Material
The code provided by this repository is part of the supplementary material for the scientific article "Impact of the used solvent on the reconstitution efficiency of evaporated biosamples for untargeted metabolomics studies" by Sascha K. Manier et al., Metabolomics, 2020. See the following sections for detailed information about it.


# Dependencies
The code deployed by this repository requires the following packages:

- XCMS (Bioconductor)
- tidyverse (CRAN)

If not installed you can use the following lines to install them:

	install.packages("tidyverse")
	
	if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("xcms")
	
# evaluate*.R
## General
These scripts were used to evaluate the peak area of the monitored compounds, as well as the total feature count. Scripts that contain an "H" in their name were used for those analyses in which HILIC was used, scripts that contain a "PH" were used for those that were obtained after usage of a PhenylHexyl column. The suffixes "pos" and "neg" describe their usage for each ionization mode.

## Usage
To use the scripts as deployed in this repository sort the files by their chromatographycal properties and then by ionisation mode. Place the matching script in the same folder were the files are grouped in folders by their reconstitution mixture and run it.

## Result
You will obtain the following graphics as PDF files:

- IS\_Area.pdf
- Creatinine\_Area.pdf (only analyses using positive ionization mode)
- GPCho\_16\_0\_16\_0\_Area.pdf (only analyses using positive ionization mode)
- Glucose\_Area.pdf (only analyses using negative ionization mode) 
- Palmitic\_Acid\_Area.pdf (only analyses using negative ionization mode)
- Peaknumbers.pdf
- Trp-d5\_EICs.pdf

The first graphics contain boxplots of the peak areas of the investigated compounds. The second graphic is a also displaying boxplots but those of the total amount of features, that were found during peak picking by XCMS. The last graphic shows every peak of >tryptophan-d<sub>5</sub> in the investigated data set. Statistics were conducted by comparing variabilities in the data set using a one-way ANOVA. Additionally Welch two sample t-test has been used to compare the means of each group against the group "QC".

# Data files
Data files can be found on MetaboLights via the study identifier [MTBLS1219](https://www.ebi.ac.uk/metabolights/MTBLS1219)
