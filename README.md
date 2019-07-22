# Code to analyze multi-lab LOD / LOQ study data

Christopher M. Merkes, Katy E. Klymus, Michael J. Allison, Caren Goldberg, Caren C. Helbing, Margaret E. Hunter, Craig A. Jackson, Richard F. Lance, Anna M. Mangan, Emy M. Monroe, Antoinette J. Piaggio, Joel P. Stokdyk, Chris C. Wilson, Catherine Richter

This code supports the manuscript "Reporting the limits of detection (LOD) and quantification (LOQ) for environmental DNA assays". The code analyzes outputs generated by another more generic code to analyze raw qPCR data for LOD and LOQ. The generic LOD / LOQ calculator code can be found at https://github.com/cmerkes/qPCR_LOD_Calc or https://doi.org/10.5066/P9GT00GB. This code requires the generic LOD / LOQ calculator script to read in and prepare the data, with additional script runs as described in comments within this code. We wrote this code to specifically analyze our data and generate the figures for the above-mentioned manuscript, however it could be adaptped by a savvy programmer to run similar analyses.

## Suggested Citation

Merkes CM, Klymus KE, Allison MJ, Goldberg C, Helbing CC, Hunter ME, Jackson CA, Lance RF, Mangan AM, Monroe EM, Piaggio AJ, Stokdyk JP, Wilson CC, Richter C. (2019) Code to analyze multi-lab LOD / LOQ study data. R Script. Available at: https://github.com/cmerkes/LOD_Analysis. DOI: https://doi.org/10.5066/P9G4MPVQ. Date Accessed: <DATE>
  
## Code files

This repository contains the following files:
- `README.md`: This file
- `LICENSE`: The standard USGS software license
- `LoD-Analysis.R`: The R code used to analyze our data and generate figures

## Data

The data used in this study are available at DOI: https://doi.org/10.5066/P9AKHU1R

## Contact for questions about the code

Primary code developer: Chris Merkes (cmerkes@usgs.gov)

## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www2.usgs.gov/visual-id/credit_usgs.html#copyright/).


This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

This software is provided "AS IS".
