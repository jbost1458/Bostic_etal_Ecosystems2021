# Bostic_etal_Ecosystems2021
Data and Scripts for Bostic et al. 2021 Ecosystems Manuscript

Data files:

isotope_NO3_Bostic2021.csv - water quality (d18O, D17O, NO3) data corresponding to all surface water samples


Column names and descriptions:

Site = Name of watershed, abbreviated. See Table 1 in Bostic et al. 2021 for full names
Date = Date of sample collection
d18O = delta 18O of NO3 (per mille)
d17O = delta 17O of NO3 (per mille)
NO3_Total = Total No3 concentrations (mg N/L)

Folder with Discharge Weighted Mean Concentrations ("conc" folder)
Columns 1:800 = 800 bootstrapped samples of water-year discharge weighted mean concentration of total NO3
Rows = Specific water year corresponding to bootstrap samples

Folder with Water Year Total Nitrate fluxes ("flux" folder)
Columns 1:800 = 800 bootstrapped samples of water-year fluxes in kg/d of total NO3
Rows = Specific water year corresponding to bootstrap samples


R scripts:
Bostic_etal_2021UncertaintyScript.R
-R script used to estimate mean and uncertainty of D17O, Flow-weighted unprocessed Atmospheric Nitrate, and Processing Efficiency 
values for Water Years 2016-2017
