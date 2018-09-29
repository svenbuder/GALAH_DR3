Repository for diagnostics of GALAH DR3
---------------------------------------

AUTHOR
-------

Sven Buder (**SB**, MPIA, WG4): buder at mpia dot de

OVERVIEW
--------

This repository includes the following subsets for GALAH DR3

**Input**:

This part is based on the input from WG3 (**sobject_iraf_53.fits**).

This file has been extended with ADQL-based x-matches with 2MASS, <i>Gaia</i> DR2 (incl. Bayesian distances from Bailer-Jones+2018), AllWISE, PanSTARRS DR1, and asteroseismic data from K2.

![alt text](input/figures/parallax_uncertainties.png "Overview of parallax quality as well as importance of use of Bayesian distances" )

** Performance for Gaia FGK Benchmark stars (Jofre+2018) and stars with asteroseismic information**:

We compare with the values from Vizier/III/281/gbs (Jofre+2018, arXiv:1808.09778) as well as the difference in the performance for stars when also nu_max values are available to constrain logg:

![alt text](gbs/figures/gbs_performance_lbol.png "Performance of the GALAH bolometric pipeline for the Gaia FGK benchmark stars")

![alt text](seis/figures/seis_comparison_3setups.png "Performence of 3 different pipelines for asteroseismic sample")
![alt text](seis/figures/seismic_sample_delta_lbol.png "Performance of bolometric pipeline for asteroseismic sample (w.r.t. asteroseismic pipeline)")
![alt text](seis/figures/seis_setup_difference_lbol.png "Histogram of differences for stellar parameters of bolometric and asteroseismic pipeline")

** Clusters**:

![alt text](clusters/figures/CMD_Kiel_FehAge.png "Overview of open and globular clusters")

**Dynamics**:

![alt text](dynamics/figures/action_overview_clean_all.png "Toomre diagram and Action plot for (clean) GALAH sample")