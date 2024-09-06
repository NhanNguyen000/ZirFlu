# CiiM | Influenza project, ZirFlu cohort

## Description
This repository for analysing the metabolomics data from ZirFlu (an influenza-vaccinated cohort) in R.

## Instructions
- The **processedData** folder stores intermediate data files which are byproducts of the data analysis and are used to produce plots.
- The **scripts** is the R script to run the data analysis. The filenames of the script files contain sequence numbers that determine their excutive order.
- The **reference** stores reference data that are used in the R scripts
- The **output** folder store the outcome of the data analysis, including:
   - Plots for main and supplement figures (.png and .svg formats)
   - Table "endogenous_metabolites.txt" contains the list of endogenous metabolites used in the metabolite analysis
   - Excel files: the outcome of linear analyses between metabolite concentrations and disease condition or HAI titers 
   - .txt files: the list of possible compound ID for statisitically significant metabolites, i.e. "meboDE", in corresponding analysis. These compound ID lists are used to run the pathways analysis in Metabo Analyst website ([metaboanalyst.ca](https://www.metaboanalyst.ca))
   - The **pathwayAnalsyses_fromMetaboAnalyst** folder store the outcome of pathway analysis in Metabo Analyst website ([metaboanalyst.ca](https://www.metaboanalyst.ca))
- The **dev** folder stores data, code, outcome at the development phase, which is a preliminary stage before the final analysis.

## Related publication
TBC. Publication DOI: To be added.
Link to our GitHub group repository: [CiiM-Bioinformatics-group/ZirFlu/](https://github.com/CiiM-Bioinformatics-group/ZirFlu/)

## Contact:
LinkedIn:	[nhannguyen](https://www.linkedin.com/in/nhannguyen1412) | ORCID: [0000-0001-8720-1195](https://orcid.org/0000-0001-8720-1195)
