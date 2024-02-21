Reproducing analysis in Reversal of Transient European Summer Precipitation Decline in a Stabilising Climate by Dittus et al., Revision 1
November 2023

----- DATA DESCRIPTION --------

Variables used: 
- tas: surface air temperature
- pr: precipitation
- huss: surface specific humidity
- msftyz: overturning streamfunction
- psl: sea-level pressure 
- ua@850 hpa: 850 hPa zonal winds
- siconca: sea-ice concentration

For peer review, data are labelled by model internal suite-id, see keys below (e.g. 'bw848', etc.)
The concatenated files are provided for peer-review (see prep1 below).

The CEDA data archive is structured as follows: 

Filenames of fixed concentration runs (for full dataset to be published on CEDA):
r1i1p1f2 (suite-id: u-bq777): 2014 in historical 
r2i1p1f2 (suite-id: u-bu607): 2040 in SSP3-7.0
r3i1p1f2 (suite-id: u-bw848): 2020 in SSP3-7.0 
r4i1p1f2 (suite-id: u-bw987): 2030 in SSP1-1.9
r5i1p1f2 (suite-id: u-bz227): 2025 in SSP3-7.0 
r6i1p1f2 (suite-id: u-cd269): 2025 in SSP2-4.5 

All simulations are branched from the first initial condition realisation of the parent ScenarioMIP experiment (r1i1p1f2 in CMIP6 nomenclature).

Filenaming convention example: 
[varname]_[miptable]_UKESM1-0-LL_futControl_r5i1p1f2_gn_245001-252412.nc
Note futControl is used to reflect the parallels to piControl (fixed external forcings, but future conditions).


----- PROCESSING & ANALYSIS --------

## Model data pre-processing and analysis: Steps prep1-3 produce the processed files needed for the plotting scripts

 * prep1 prep_intermediate_files.py: produce concatenated files in a consistent format, extract specified pressure level (concatenated files are uploaded for peer-review, but will not be uploaded in final repository - concatenation can be performed on model archived files from CEDA)

 * prep2 Jet_lat.ipynb: calculate jet time series and save to netcdf

 * prep3 Regional_preproc.ipynb: calculate precipitation for different seasons and IPCC regions following Iturbide et al., save as timeseries input for later scripts 

## Figures: 

 * Figure 1:
	- This figure can be reproduced by running the script named: Figure1.ipynb	

 * Figure 2:
    - This figure can be reproduced by running the script named: Figure2.ipynb	

 * Figures 3 and 4: 

	- These figures can be reproduced by running the script named: Figures3and4.ipynb
	- Input:  the regional time series and jet time series, as well as the intermediate model data structure / files (output of pp1, pp2 and pp3)

 * Supplementary figures:
   - Figure S1: Figures3and4.ipynb using global domain
   - Figure S2: Included in Figure1.ipynb
   - Figure S3: SAT_AMOC_regression.ipynb
   - Figure S4: As Figure 2 with different trend length
   - Figure S5 & S6: Figures3and4.ipynb with different values for length and season

