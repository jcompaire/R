## An R repository to perform a time series modeling based on the Autoregressive Integrated Moving Average (ARIMA) models

This repository contains R scripts to replicate the analyses performed in the manuscript 
"***Time series modeling of coastal fishery landings on the Southwestern Atlantic shelf***", 
which has been submitted to [Fishery Bulletin](https://spo.nmfs.noaa.gov/fb.htm). The code can be used not only to compute statistical analyses but also to generate all figures and outputs presented in the manuscript. Our goal is to increase transparency and reproducibility in research.

### Overwiew
We investigated the dynamics of fishery ladings for the three fish species most caught at the Rio de la Plata estuary and its maritime front (Southwestern Atlantic Ocean) using Autoregressive Integrated Moving Average (ARIMA) models. Our results prove that the selected ARIMA models are able to capture the dynamics of fishery landings and provide a short-term reliable prediction of catches for the three species. 

### Authors
This work was developed through a collaboration between [Jesus C. Compaire](https://www.researchgate.net/profile/Jesus-Compaire), [E. Marcelo Acha](https://www.researchgate.net/profile/Marcelo-Acha), [Diego Moreira](https://www.researchgate.net/profile/Diego-Moreira-3) and [Claudia G. Simionato](https://www.researchgate.net/profile/Claudia-Simionato). The R code was developed and is maintained by Jesus C. Compaire.

### Files description and workflow
- ***time_series_aucfz_data.RData*** contains the time series of monthly commercial landings for the 1996-2021 period for Cynoscion guatucupa, Micropogonias furnieri and Merluccius hubbsi.
- ***time_series_aucfz_functions.R*** contains the algorithms and functions necessary to run the file "time_series_aucfz_code.R" .
- ***time_series_aucfz_code.R*** contains the code to load the dataset, perform the statistical analyses and get the outputs.

The workflow to replicate the analyses is pretty easy: just ensure that you have in the same path the three previously mentioned files prior to running the last one. Note: necessary packages are listed at the beginning of thes file *time_series_aufc_code.R*.

## References

Compaire, J.C., Acha, E.M., Moreira, D. & Simionato C.G. (2023). Time series modeling of coastal fishery landings on the Southwestern Atlantic shelf (submitted to Fishery Bulletin)
