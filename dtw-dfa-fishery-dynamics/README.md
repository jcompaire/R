Multiscale Environmental Forcing of Multispecies Fishery Dynamics

This repository contains the data and code to replicate the statistical analyses and figures presented in the manuscript:

    Compaire, J.C., Irigoyen, A.J., Simionato, C.G. & Acha, E.M. (2025). Multiscale environmental forcing of multispecies fishery dynamics across climate timescales. Submitted to Fisheries Oceanography.

🌐 Interactive Project Website

The full analysis, including code execution and interactive visualizations, is available as a Quarto website:
[https://github.com/jcompaire/R/dtw-dfa-fishery-dynamics]
🔬 Overview

This research explores the multiscale environmental drivers affecting fishery dynamics in the Argentine-Uruguayan Common Fishing Zone (AUCFZ). First, we applied a time-series alignment approach to account for potential temporal asynchronies and a latent factor modeling framework to identify common patterns across the multivariate landing series. Then,  resulting patterns were used as response variables in lagged generalized linear and additive models, incorporating a comprehensive suite of environmental predictors as covariates:

    Dynamic Time Warping (DTW) & Dynamic Factor Analysis (DFA): To identify common patterns and trends in multispecies catch data across different timescales.

    Generalized Linear Models (GLM) & Generalized Additive Models (GAM): To evaluate the relationship between environmental variables and fishery variability.

📁 Repository Structure

    1_dtw_dfa_analyses.qmd: Detailed workflow for DTW and DFA analyses.

    2_glm_gam_analyses.qmd: Detailed workflow for GLM and GAM analyses.
    
    dtw_dfa_glm_gam/dtw_dfa_glm_gam_aucfz_functions.R: Auxiliary R functions used throughout the analysis.

    datasets/: Contains the necessary .nc, .tiff, and .RData files required to reproduce the results.

    docs/: The rendered HTML files for the project website.

🛠️ Reproducibility

To replicate this research locally, you will need R installed along with the Quarto CLI.

    Clone the repository: git clone [https://github.com/jcompaire/R/dtw-dfa-fishery-dynamics.git](https://github.com/jcompaire/R/dtw-dfa-fishery-dynamics.git)

    Open the project in RStudio or your preferred IDE.

    Render the entire project using the terminal: quarto render.

👤 Authors

    [Jesus C. Compaire](https://www.researchgate.net/profile/Jesus-Compaire) | Lead Developer: wrote and maintains the repository and its contents.
    Alejo J. Irigoyen
    Claudia G. Simionato
    Eduardo M. Acha

This work is part of our commitment to Open Science. All code is shared under a MIT License.

