# Supplementary Material
This Respository contains supplementary data for a submitted manuskript. Exact title and more information will be provided after acceptance.

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
These scripts were used to evaluate the peak area of the internal standard, as well as the total feature count. Scripts that contain an "H" in their name were used for those analyses in which HILIC was used, scripts that contain a "PH" were used for those that were obtained after usage of a PhenylHexyl column. The suffixes "pos" and "neg" describe their usage for each ionization mode.

## Usage
To use the scripts as deployed in this repository sort the files by their chromatographycal properties and then by ionisation mode. Place the matching script in the same folder were the files are grouped in folders by their reconstitution mixture and run it.

## Result
You will obtain three graphics as PDF files:

- IS_Area.pdf
- Peaknumbers.pdf
- Trp-d5_EICs.pdf

The first graphic contains boxplots of the peak areas of the investigated internal standard Tryptophan-d<sub>5</sub>. The second graphic is a barplot displaying the total amount of features, that was found during peak picking by XCMS. The last graphic shows every peak of >tryptophan-d<sub>5</sub> in the investigated data set.

# Data files
Data files can be found on MetaboLights via the study identifier [MTBLS1219](https://www.ebi.ac.uk/metabolights/MTBLS1219)
