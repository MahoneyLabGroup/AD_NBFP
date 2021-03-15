# Alzheimer's Disease Network-Based Gene Prioritization

This repository contains all the code necessary to replicate the analysis performed in the manuscript ["System-level analysis of Alzheimer's disease prioritizes candidate genes for neurodegeneration"](https://www.frontiersin.org/articles/10.3389/fgene.2021.625246/abstract). We have provided gene-level summary statistics for the large MetaGWAS and ADNI datasets. To get access to SNP-level summary statistics for the MetaGWAS, please visit the website for the [Complex Trait Genetics lab](AD_sumstats_Jansenetal_2019sept.txt.gz). Access to ADNI data is limited to those who apply and receive access to the data. Please visit the [ADNI](http://adni.loni.usc.edu/data-samples/access-data/) website for more information on how to apply.

# Organization
---
* The methods.Rmd file contains all the code to replicate the analysis. Please read the text as it contains information on files that are available to be read in as well as how to run the SVMs on your own.
* The *code* directory contains all the custom source code needed to run the analysis. 
* The *data* directory contains the gene-level summary statistics for the ADNI and MetaGWAS data. It is also the location where the genetic networks will be stored once downloaded from [HumanBase](https://hb.flatironinstitute.org).
* The *mwas* directory will contain the results from the SVMs once that code is run. The full results are too large to upload to GitHub, though we have included a final scored and merged results table for both brain tissues which can be found in the *results* directory along with results from REVIGO that are used to construct the pie charts.

# Required Packages
---
The latest version of this workflow was tested in R/4.0.3. 

An internet is required for this workflow as it uses Biomart and `gprofiler2` to retrieve information about genes. The full list of required packages can be installed using the first chunk in the `methods.Rmd` file.
