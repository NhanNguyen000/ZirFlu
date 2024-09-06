# CiiM | Influenza project, ZirFlu cohort

## Description
This repository for analysing the metabolomics data from ZirFlu (an influenza-vaccinated cohort) in R.

## Instructions
- The **processedData** folder stores intermediate data files which are byproducts of the data analysis and are used to produce plots.
- The **scripts** is the R script to run the data analysis. The filenames of the script files contain sequence numbers that determine their excutive order.
- The **reference** stores reference data that are used in the R scripts
- The **output** folder store the outcome of the data analysis, including:
   - (1) plots (.png and .svg formats) for main and supplement figures; 
   - (2) "endogenous_metabolites.txt" table for the list of endogenous metabolites used in the metabolite analysis; 
   - (3) the excel files for the outcome of linear analyses between metaboliti concentrations and disease condition or HAI titers;  
   - (4) all .txt files start with "meboDE" are the list of possible compound ID for statisitically significant metabolites in corresponding analysis. These compound ID lists are used to run the enrichment analysis in Metabo Analyst website (metaboanalyst.ca)
- The **dev** folder stores data, code, outcome at the development phase, which is a preliminary stage before the final analysis.

## Related publication
TBC. Publication DOI: To be added.

## Contact:
LinkedIn:	[nhannguyen](https://www.linkedin.com/in/nhannguyen1412) | ORCID: [0000-0001-8720-1195](https://orcid.org/0000-0001-8720-1195)
