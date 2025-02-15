# ibex_gloworm_model

Eleanor Dickinson November 2023

This repository contains data and code to reproduce the results presented in a paper by Dickinson et al. 2024 published in Royal Society Open Science, titled “Host movement dominates the predicted effects of climate change on parasite transmission between wild and domestic mountain ungulates”. https://doi.org/10.1098/rsos.230469. 

The following folders contain: 
1. host_elevation

     1.1 hostelevation.R: R script for predicting host elevation.

   1.2 avg_elev_week_2017-19.csv: raw data for host elevation.

   1.3 RCP85_2065-2095_365.csv: Projected climate data from the high emissions scenario (Representative Concentration Pathway 8.5; RCP 8.5) from the HADGEM-ES model output for a 30-year period (2066 - 2095) (Martin et al., 2011).
2. parameters

    2.1  belledonnedata_FULL.csv: the full data for the study area included in the study.
   
    2.2  Levionazages.csv: age data for Alpine ibex (Brivio et al., 2019).

   2.3 Gruner_validation: larval pasture counts reported in Gruner et al., 2008. 
3. code

     3.1 glowormfunctions_dickinson.R: R script file containing the functions used in the model used to simulate parasite dynamics.
   
     3.2 glowormrun_dickinson.R: R script containing the code to run the model simulations.

     3.3 glowormoutput_dicinson.R: R script containing code to summarize the output, calculate host exposure, and run the validation.

   3.4 glowormclimate_dickinson.R: R script containing the code for the historic and projected climate scenarios. 

     




References

Brivio, F., Zurmühl, M., Grignolio, S., von Hardenberg, A., Apollonio, M., & Ciuti, S. (2019) Forecasting the response to global warming in a heat-sensitive species. Scientific Reports, 9,  3048. 

Gruner, L., Sauvé, C., Boulard, C., & Calamel, M. (2006). Analysis of the relationship between land use and the parasitism of sheep during their transhumance. Animal Research, 55 (3), 177–188.

Martin, G.M., Bellouin, N., Collins, W.J., Culverwell, I.D., Halloran, P.R., Hardiman, S.C., et al. (2011) The HadGEM2 family of Met Office Unified Model climate configurations. Geoscientific Model Development, 4 (3), 723–757. 
