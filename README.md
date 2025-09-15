# Spatially Explicit Capture--Recapture Simulations

[![DOI](https://zenodo.org/badge/785481834.svg)](https://doi.org/10.5281/zenodo.17127878)

This site holds R code and results for simulations to assess the robustness of spatially explicit capture--recapture estimates of population density to various breaches of model assumptions. Results are used in [Efford (2025)](https://murrayefford.github.io/SECRbook/).

All simulations were conducted in the R package **secrdesign** and SECR models were fitted by maximum likelihood in R package **secr**. Some shared code and other detail is given in the file 'setup.R' in the root folder.

Simulations used the New Zealand eScience Infrastructure (NeSI) high performance computing facilities URL https://www.nesi.org.nz. Prior to 2025, batch operation on that platform did not allow R graphics, so final rendering of the .rmd files was performed on a desktop machine using stored simulation results (the .RDS files in each folder).

Follow the HTML link for results and the RMD link for R code.

| HTML | RMD | Simulation| 
|:-|:--|:--------------------|
| [ARR]  | [ARR rmd]  | Size of array |
| [CLO]  | [CLO rmd]  | Failure of closure due to mortality and recruitment |
| [DNC]  | [DNC rmd]  | Failure of closure due to dispersal |
| [LDF]  | [LDF rmd]  | Latent detection field (Stevenson et al. 2021) |
| [MID]  | [MID rmd]  | Mis-identification |
| [Mb]   | [Mb rmd]   | Behavioural effects |
| [Mh]   | [Mh rmd]   | Individual heterogeneity |
| [MhS]   | [MhS rmd]   | Spatial trend  in presence of individual heterogeneity |
| [MhT]   | [MhT rmd]   | Temporal trend in presence of individual heterogeneity |
| [Mt]   | [Mt rmd]   | Temporal variation |
| [OU]   | [OU rmd]   | Serial correlation of location within home range |
| [RSF]  | [RSF rmd]  | Habitat-driven variation in related detection<sup>1</sup> |
| [SARE] | [SARE rmd] | Spatially autocorrelated random effects<sup>1</sup> |
| [STR]  | [STR rmd]  | Stratification and spatial variation in density, effort and detection |

1. RSF and SARE are different aspects of the same phenomenon; RSF relates to Royle et al. (2013) and Efford (2014), whereas SARE relates to Moqanaki et al. (2021) and Dey et al. (2023).

Please report Issues on this GitHub page. 

### References

Dey, S., Moqanaki, E., Milleret, C., Dupont, P., Tourani, M. and Bischof, R. 2023. Modelling Spatially autocorrelated detection probabilities in spatial capture--recapture using random effects. Ecological Modelling 479: 110324.

Efford, M. G. 2014. Bias from heterogeneous usage of space in spatially explicit capture--recapture analyses. Methods in Ecology and Evolution 5: 599--602.

Moqanaki, E. M., Milleret, C. Tourani, M., Dupont, P. and Bischof, R. 2021. Consequences of ignoring variable and spatially autocorrelated detection probability in spatial capture--recapture. Landscape Ecology 36: 2879--95. 

Royle, J. A., Chandler, R.B., Sun, C.C., and Fuller, A. K. 2013. Integrating resource selection information with spatial capture--recapture. Methods in Ecology and Evolution 4: 520--30.

Stevenson, B.C., Fewster, R. M. and Sharma, K. 2021. Spatial correlation structures for detections of individuals in spatial capture--recapture models. Biometrics 78: 963--73.

Last updated: 16 September 2025

[ARR]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/ARR/secr-simulations-ARR.html
[CLO]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/CLO/secr-simulations-CLO.html
[DNC]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/DNC/secr-simulations-DNC.html
[LDF]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/LDF/secr-simulations-LDF.html
[MID]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/MID/secr-simulations-MID.html
[Mb]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/Mb/secr-simulations-Mb.html
[Mh]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/Mh/secr-simulations-Mh.html
[MhS]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/MhS/secr-simulations-MhS.html
[MhT]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/MhT/secr-simulations-MhT.html
[Mt]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/Mt/secr-simulations-Mt.html
[OU]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/OU/secr-simulations-OU.html
[RSF]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/RSF/secr-simulations-RSF.html
[SARE]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/SARE/secr-simulations-SARE.html
[STR]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/STR/secr-simulations-STR.html

[ARR rmd]: SARE/secr-simulations-ARR.rmd
[CLO rmd]: CLO/secr-simulations-CLO.rmd
[DNC rmd]: DNC/secr-simulations-DNC.rmd
[LDF rmd]: LDF/secr-simulations-LDF.rmd
[MID rmd]: DNC/secr-simulations-MID.rmd
[Mb rmd]: Mb/secr-simulations-Mb.rmd
[Mh rmd]: Mh/secr-simulations-Mh.rmd
[MhS rmd]: MhS/secr-simulations-MhS.rmd
[MhT rmd]: MhT/secr-simulations-MhT.rmd
[Mt rmd]: Mt/secr-simulations-Mt.rmd
[OU rmd]: OU/secr-simulations-OU.rmd
[RSF rmd]: RSF/secr-simulations-RSF.rmd
[SARE rmd]: SARE/secr-simulations-SARE.rmd
[STR rmd]: STR/secr-simulations-STR.rmd
