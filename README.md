# gpu-SAC
GPU-SAC is a GPU implementation of the Sacramento Soil Moisture Accounting (SAC-SMA) model. This hydrological model takes into account soil moisture content and forecast precipitation as well as temperatures and it simulates runoff.

Schematic of the model is found in the following image ([UC Irvine](http://chrs.web.uci.edu/research/hydrologic_predictions/activities07.html))

![Image of Yaktocat](http://chrs.web.uci.edu/images/sac_smal_medium.jpg)

More information can be found [here](http://www.cbrfc.noaa.gov/wsup/sac_sm/cbrfc_sacsma_101_20140731.pdf) and [here](http://www.manureadvisorysystem.wi.gov/app/SACmodel).

## Files
- __GPU-SAC.c__ Testing application. It can be compiled using the Makefile
- __LICENSE.txt__ License of the project
- __Makefile.mak__ Contains the rules for the __make__ command. 
- __src__
 - __sac_define.h__ General definitions and configuration
 - __sac_structs.h__ Data-types and structures
 - __sac_engine.*__ Main functions for the SAC-SMA model (CPU)
 - __sac_g_engine.*__ Main functions for the SAC-SMA model (GPU)
 - __sac_io.h__ input/output functions
 - __sac_utilities.h__ Miscelanea
 - 
# Project Website
Click [here](http://alacasta.github.io/gpu-SAC/) for the GPU-SAC project website. 

## License
GPU-SAC is an OSS licensed under [MIT License](https://github.com/alacasta/gpu-SAC/blob/master/LICENSE.txt).





