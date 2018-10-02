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

**Performance for Gaia FGK Benchmark stars (Jofre+2018) and stars with asteroseismic information**:

We compare with the values from Vizier/III/281/gbs (Jofre+2018, arXiv:1808.09778) as well as the difference in the performance for stars when also nu_max values are available to constrain logg:

We do not see any significant biases for Teff and logg, but had to correct and [Fe/H] bias of -0.1 (underestimated [Fe/H]).

![alt text](gbs/figures/gbs_performance_lbol.png "Performance of the GALAH bolometric pipeline for the Gaia FGK benchmark stars")

The overlap of the stars with asteroseismic parameters and parallaxes show that the bolometric pipeline (middle panels) performs significantly better than the pipeline without additional non-spectroscopic information (left panels). Especially the Red Clump stars show an outstanding agreement with the pipeline results which also used nu_max values.

![alt text](seis/figures/seis_comparison_3setups.png "Performence of 3 different pipelines for asteroseismic sample")

The scatter in the differences between asteroseismic and bolometric pipeline for logg is mainly driven by the red clump stars, which have to be further investigated, but are significantly better than for GALAH DR2.

![alt text](seis/figures/seismic_sample_delta_lbol.png "Performance of bolometric pipeline for asteroseismic sample (w.r.t. asteroseismic pipeline)")

The biases between asteroseismic and bolometric pipelines for the stellar parameters are not significant but show a trend of lower Teff (-21K), lower logg (-0.03dex), and lower [Fe/H] (-0.03dex) than the values from the pipeline using asteroseismic values/relations for logg.

![alt text](seis/figures/seis_setup_difference_lbol.png "Histogram of differences for stellar parameters of bolometric and asteroseismic pipeline")

**Abundance zeropoints**

The current abundance zeropoints for the Sun and Arcturus are very small for light elements but need further analyses for the heavy elements.

![alt text](abundance_zeropoints/figures/abundance_zeropoints.png "Abundance zeropoints for the Sun and Arcturus")

**Clusters**:

Plotted are all globular and open clusters members (based on the membership analysis by Janez Kos from 180507 and if SME parameters available):

![alt text](clusters/figures/CMD_Kiel_FehAge.png "Overview of open and globular clusters")

![alt text](clusters/figures/CMD_Kiel_FehAge_M67.png "Overview for M 67")

Left: CMD, middle: Kiel diagram, right: 2D-histogram of ages and [Fe/H] of all cluster stars

**Dynamics**:

![alt text](dynamics/figures/action_overview_clean_all.png "Toomre diagram and Action plot for (clean) GALAH sample")