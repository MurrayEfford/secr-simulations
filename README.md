# Spatially Explicit Capture-Recapture Simulations

This site holds R code and results for simulations to assess the robustness of spatially explicit capture--recapture estimates of population density to various breaches of model assumptions.

All simulations were conducted in the R package **secrdesign** and SECR models were fitted by maximum likelihood in R package **secr**.

Follow the HTML link for results and the RMD link for R code.

| HTML | RMD | Simulation| 
|:-|:--|:--------------------|
| [CLO]  | [CLO rmd]  | Failure of closure due to mortality and recruitment |
| [DNC]  | [DNC rmd]  | Failure of closure due to dispersal |
| [MID]  | [MID rmd]  | Misidentification |
| [Mt]   | [Mt rmd]   | Temporal variation |
| [Mb]   | [Mb rmd]   | Behavioural effects |
| [Mh]   | [Mh rmd]   | Individual heterogeneity |
| [OU]   | [OU rmd]   | Serial correlation of location within home range |
| [SARE] | [SARE rmd] | Spatially autocorrelated random effects^1^ |
| [LDF]  | [LDF rmd]  | Latent detection field Stevenson et al. (2021) |

1. SARE includes the scenarios considered by Royle et al. (2013), Efford (2014), Moqanaki et al. (2021), and Dey et al. (2023)

Please report Issues on this GitHub page. 

### References

Dey et al. 2023  
Efford 2014
Moqanaki et al. 2021  
Royle et al. 2013
Stevenson et al. 2021  


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
