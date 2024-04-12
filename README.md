# Spatially Explicit Capture-Recapture Simulations

SITE UNDER DEVELOPMENT - RESULTS WILL CHANGE

This site holds R code and results for simulations to assess the robustness of spatially explicit capture--recapture estimates of population density to various breaches of model assumptions.

All simulations were conducted in the R package **secrdesign** and SECR models were fitted by maximum likelihood in R package **secr**.

Follow the HTML link for results and the RMD link for R code.

| HTML | RMD | Simulation| 
|:-|:--|:--------------------|
| [CLO]  | [CLO rmd]  | Failure of closure due to mortality and recruitment |
| [DNC]  | [DNC rmd]  | Failure of closure due to dispersal |
| [MID]  | [MID rmd]  | Mis-identification |
| [Mt]   | [Mt rmd]   | Temporal variation |
| [Mb]   | [Mb rmd]   | Behavioural effects |
| [Mh]   | [Mh rmd]   | Individual heterogeneity |
| [OU]   | [OU rmd]   | Serial correlation of location within home range |
| [SARE] | [SARE rmd] | Spatially autocorrelated random effects^1^ |
| [LDF]  | [LDF rmd]  | Latent detection field Stevenson et al. (2021) |

1. SARE includes the scenarios considered by Royle et al. (2013), Efford (2014), Moqanaki et al. (2021), and Dey et al. (2023)

Please report Issues on this GitHub page. 

### References

Dey, S., E. Moqanaki, C. Milleret, P. Dupont, M. Tourani, and R. Bischof. 2023. Modelling Spatially Autocorrelated Detection Probabilities in Spatial Capture-Recapture Using Random Effects. Ecological Modelling 479: 110324.

Efford, M. G. 2014. Bias from Heterogeneous Usage of Space in Spatially Explicit Capture-Recapture Analyses. Methods in Ecology and Evolution 5: 599--602.

Moqanaki, E. M., Milleret, C. Tourani, M., Dupont, P. and Bischof, R. 2021. Consequences of Ignoring Variable and Spatially Autocorrelated Detection Probability in Spatial Capture-Recapture. Landscape Ecology 36: 2879--95. 

Royle, J. A., Chandler, R.B., Sun, C.C., and Fuller, A. K. 2013. Integrating Resource Selection Information with Spatial Capture–Recapture.” Methods in Ecology and Evolution 4: 520--30.

Stevenson, B.C., Fewster, R. M. and Sharma, K. 2021. Spatial Correlation Structures for Detections of Individuals in Spatial Capture–Recapture Models. Biometrics 78: 963--73.

[CLO]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/CLO/secr-simulations-CLO.html
[DNC]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/DNC/secr-simulations-DNC.html
[MID]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/MID/secr-simulations-MID.html
[LDF]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/LDF/secr-simulations-LDF.html
[Mb]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/Mb/secr-simulations-Mb.html
[Mt]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/Mt/secr-simulations-Mt.html
[Mh]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/Mh/secr-simulations-Mh.html
[OU]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/OU/secr-simulations-OU.html
[SARE]: https://htmlpreview.github.io/?https://github.com/MurrayEfford/secr-simulations/blob/main/SARE/secr-simulations-SARE.html

[CLO rmd]: CLO/secr-simulations-CLO.rmd
[DNC rmd]: DNC/secr-simulations-DNC.rmd
[MID rmd]: DNC/secr-simulations-MID.rmd
[LDF rmd]: LDF/secr-simulations-LDF.rmd
[Mb rmd]: Mb/secr-simulations-Mb.rmd
[Mt rmd]: Mt/secr-simulations-Mt.rmd
[Mh rmd]: Mh/secr-simulations-Mh.rmd
[OU rmd]: OU/secr-simulations-OU.rmd
[SARE rmd]: SARE/secr-simulations-SARE.rmd
