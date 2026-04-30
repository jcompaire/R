## An R repository for running Generalized Linear Models (GLMs) and Generalized Additive Models (GAMs)

This repository contains R scripts to replicate the analyses performed in the manuscript 
[***Modeling environmental effects on fishery landings: A case study of whitemouth croaker (_Micropogonias furnieri_) in the Southwestern
Atlantic shelf***](https://doi.org/10.1016/j.ecss.2024.108806). The code can be used not only to compute statistical analyses but also to generate all figures and outputs presented in the manuscript. Our goal is to increase transparency and reproducibility in research.

### Overwiew
We investigated the effect of environmental variability on commercial landings of the whitemouth croaker (_Micropogonias furnieri_) caught at the Rio de la Plata estuary and its maritime front (Southwestern Atlantic Ocean) using Generalized Linear Models (GLMs) and Generalized Additive Models (GAMs). Our results prove that GAMs outperform GLMs in capturing the underlying dynamics when assessing the impact of environmental variability on fisheries.

### Authors
This work was conducted by [Jesus C. Compaire](https://www.researchgate.net/profile/Jesus-Compaire), [Claudia G. Simionato](https://www.researchgate.net/profile/Claudia-Simionato), [Diego Moreira](https://www.researchgate.net/profile/Diego-Moreira-3) and [E. Marcelo Acha](https://www.researchgate.net/profile/Marcelo-Acha). Jesus C. Compaire wrote and updates the R scripts.

### Files description and workflow
- `elevation_AUCFZ_ETOPO2022.tiff` contains elevation data for the study area retrieved from [ETOPO](https://doi.org/10.25921/fd45-gt74) 
- `glm_gam_aucfz_data.RData` contains the time series of monthly commercial landings for the 1996-2021 period for [*Micropogonias furnieri*](https://www.fishbase.se/summary/Micropogonias-furnieri.html).
- `glm_gam_aucfz_functions.R` contains the algorithms and functions necessary to run the file "glm_gam_aucfz_code.R" .
- `glm_gam_aucfz_code.R` contains the code to load the dataset, perform the statistical analyses and get the outputs.

The workflow to replicate the analyses is pretty easy: just ensure that you have in the same path the three previously mentioned files prior to running the last one. Note: necessary packages are listed at the beginning of the file *glm_gam_aucfz_code.R*.

## References

[Compaire et al. (2024)](https://doi.org/10.1016/j.ecss.2024.108806)

```
Compaire, J.C., Simionato C.G.  Moreira, D. & Acha, E.M. (2024).
Modeling environmental effects on fishery landings: A case study of whitemouth croaker (Micropogonias furnieri)
in the Southwestern Atlantic shelf.
Estuarine, Coastal and Shelf Science, 305, 108806, https://doi.org/10.1016/j.ecss.2024.108806
```
[![DOI:10.1016/j.ecss.2024.108806](http://img.shields.io/badge/DOI-10.1016/j.ecss.2024.108806-b45f06.svg)](https://doi.org/10.1016/j.ecss.2024.108806)
