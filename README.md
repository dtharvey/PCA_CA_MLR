# PCA_CA_MLR

Materials for Pittcon 2018 talk on principal component analysis, cluster analysis, and multiple linear regression in the undergraduate analytical chemistry curriculum (part of the symposium "Data Science in the Chemistry Curriculum"). 

**Title**: Using R to Introduce Students to Principal Component Analysis, Cluster Analysis, and Multiple Linear Regression

**Abstract**: A common experiment in many undergraduate courses in analytical chemistry is the quantitative analysis of a two-component mixture by UV/Visible spectrophotometry. The analysis typically involves measuring a sample's absorbance at two wavelengths, determining the molar absorptivity for each analyte and each wavelength, and using Beer's law to solve for the concentrations of the two analytes simultaneously. Rarely is this analysis extended to more than two components or to data collected at more wavelengths than there are analytes. In this presentation, we will consider how to use the analysis of multicomponent mixtures to introduce undergraduate students to principal component analysis, cluster analysis, and multiple linear regression using the statistical programming language R.

The file slideScripts.R provides the code used to generate the figures included in the slide deck. Note that the script makes use of the readr, scatterplot3d and chemCal packages.
