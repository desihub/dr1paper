# DESI Data Release 1 

## Zenodo Data

Data files associated with the figures in DESI Collaboration et al
(2025, AJ, submitted, arXiv:2503.14745), "Data Release 1 of the Dark
Energy Spectroscopic Instrument".

The figures in this paper were generated using code and notebooks in
https://github.com/desihub/dr1paper. For convenience and for all but
the first two figures in the paper (see below), we also provide the
processed data points used to generate these figures.

## Figure 1:

This figure was generated using the full (publicly released) DR1
large-scale structure catalogs.

**Code**: `nb/DESI_Y1_butterfly.ipynb`


## Figure 2:

This figure was generated using the full (publicly released) DR1 Milky
Way Survey catalogs.

**Code**: `nb/DESI_Y1_MW.ipynb`


## Figure 3: 

**Data files**:
  * `zenodo/bright_completeness.csv` with columns `healpix_number,completeness`
  * `zenodo/dark_completeness.csv` with columns `healpix_number,completeness`
  * `zenodo/backup_completeness.csv` with columns `healpix_number,completeness`

**Code**: `bin/progressplot`


## Figure 4: number of unique tiles per night

**Data file**: `zenodo/tiles_per_night.csv` with columns
  `night,sv1,sv2,sv3,special,main_bright,main_dark,main_backup`

**Code**: `bin/count_obs`


## Figure 5: 

**Data files**:
  * `zenodo/bright_density.csv` with columns `healpix_number,density`
  * `zenodo/dark_density.csv` with columns `healpix_number,density`
  * `zenodo/backup_density.csv` with columns `healpix_number,density`
  * `zenodo/sv_density.csv` with columns `healpix_number,density`

**Code**: `bin/plot_density_skymap`


## Figure 6: 

**Data files**:
  * `zenodo/nofz_GALAXY.csv` with columns `ZBIN,ALL,ELG,LRG,BGS`
  * `zenodo/nofz_QSO.csv` with columns `ZBIN,ALL,QSO`
  * `zenodo/nofv_STAR.csv` with columns `VBIN,ALL,MWS`

**Code**: `nb/count-dr1-targets.ipynb`

## Figure 7: 

This figure was generated using Astro Data Lab tools (see the notebook
for details).

**Code**: `nb/desi_dr1_stacked_spectra_fig7.ipynb`
