# BatsOnTheMove
Scripts and datasets used in Moshier et al. 2026 (update with link) for cleaning GBIF data, modeling species distributions, and performing niche overlap calculations.

Included scripts:
1) SDM - code used to generate SDMs using BioClim, MaxEnt, GLM, Random Forest, and weighted ensemble models
2) SDM_results - processing model results and feature importance
3) centroids - file management, used to prepare files and folders to create clipping extents and calculate centroids
4) niche_overlap - used to calculate niche overlap metrics for each species using ENMTools (Warren et al. 2021), terra (Hijmans et al. 2026), and geosphere (Hijmans et al. 2024)
5) range_characteristics - used to calculate range characteristics (total area, latitudinal extent, midpoint latitude, midpoint longitude) using sf/sp/lwsgeom (Pebesma et al. 2026 papers)
6) figure_plots - code used to generate figures in the paper and supplemental information (necessary data files also included)

Included Supplemental Datasets:
1) Supplemental Table 1: literature review information for currently recognized bat species
2) Supplemental Table 2: calculated niche overlap values (Schoener's D, Hellinger's I, and rank correlation (rho)) for each bat species included in our analysis
