# Shear-stress-project
Authors: Ruben Methorst, Sander van der Laan

This repository holds all scripts used to analyze data and create figures presented in "Exploring the Causal Inference of Shear Stress Associated DNA Methylation on Cardiovascular Risk".

# **Abstract**

**Background and aims:** Atherosclerosis is a lipid-driven inflammatory disease presumably initiated by endothelial activation. Low vascular shear stress is known for its ability to activate endothelial cells. Differential DNA methylation (DNAm) is a relatively unexplored player in atherosclerotic disease development and endothelial dysfunction. Literature search revealed that expression of 11 genes have been found to be associated with differential DNAm due to low shear stress in endothelial cells. We hypothesized a causal relationship between DNAm of shear stress associated genes in human carotid plaque and increased risk of cardiovascular disease.

**Methods:** Using Mendelian randomisation (MR) analysis, we explored the potential causal role of DNAm of shear stress associated genes on cardiovascular disease risk. We used genetic and DNAm data of 442 carotid endarterectomy derived advanced plaques from the Athero-Express Biobank Study for quantitative trait loci (QTL) discovery and performed MR analysis using these QTLs and GWAS summary statistics of coronary artery disease (CAD) and ischemic stroke (IS).
**Results:** We discovered 9 methylation QTLs in plaque for shear stress associated differential DNAm. We found no significant effect of shear stress gene promotor methylation and increased risk of CAD and IS. 


#**Scripts**

* `CpG_Filtering.R`: used to filter CpG sites for only promoter CpG sites

* `MR IS and CDH.R`: Script to perform Mendelian randomisation on our data for ischemic stroke and coronary heart disease

* `RAP function.R`: script to make regional association plot

* `ROI script.R`: script to help build your regions of interest for mQTL analysis (QTLToolKit)


--------------------------
## **The MIT License (MIT)**

### Copyright (c) 1979-2020

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
