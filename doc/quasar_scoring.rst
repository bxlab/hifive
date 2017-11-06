.. _quasar:

*****************************
QuASAR Hi-C Scoring
*****************************

HiFive's QuASAR transformation and scoring (Quality Assessment of Spatial Arrangement Reproducibility) uses the consistency between the raw interaction matrix and the correlation matrix of distance-corrected signal to determine sample quality. For all interactions less than or equal to 100 times the resolution of size, the correlation value is calculated. In order to determine the quality score for a given chromosome, the weighted average of the correlation values (weighted by the raw interaction signal) minus the unweighted correlation signal average is calculated. A genome-wide value is calculated by summing all numerators and denominators across all chromosomes prior to dividing and subtracting score components. The replicate score is calculated by finding the correlation between the weighted correlation matrices for two samples.
