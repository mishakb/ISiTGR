[![N|Solid](https://www.utdallas.edu/~jdossett/images/banner_bkgd_isitgr3.jpg)](https://www.utdallas.edu/~jdossett/images/banner_bkgd_isitgr3.jpg)

# ISiTGR Version 3.1 released in February 2020 (with python wrapper)
We introduce a new version of **I**ntegrated **S**oftware **i**n **T**esting **G**eneral **R**elativity (ISiTGR) which is a patch to the software CAMB and CosmoMC. ISiTGR is intended to test deviations from GR at cosmological scales using available cosmological data sets. While doing so, it allows for various extensions to the standard flat LCDM model. In this new release, we have a combined support for the following:  

* Dynamical dark energy parametrizations with a constant or time-dependent equation of state; 

* A consistent implementation of anisotropic shear to model massive neutrinos throughout the full formalism;

* Multiple commonly used parameterizations of modified growth (MG) parameters; 

* Functional, binned and hybrid time- and scale-dependencies for all MG parameters (expanded from previous version); 

* Spatially flat or curved backgrounds (present in previous version as well). 

Additionally, **this new version of ISiTGR provides a python wrapper**. The python wrapper extends the original CAMB's python wrapper to work with the different MG parameterizations provided in ISiTGR, allowing the user to obtained power spectra and transfer functions in the context of phenomenological modified gravity.

The description of the formalism and its implementation in the CMB code, the Integrated Sachs-Wolfe (ISW) effect, and the 3x2 point statistics as well as examples of application to current data sets, can be found in the latest paper on the arXiv. A more technical description of the implementation can be found in the documentation provided in this repository. 

## Installation (to run MCMC)
The corresponding version of ISiTGR was built based on the CosmoMC July 2019 version. If you install ISiTGR from the principal folder, you will be allowed to do MCMC sampling to constraint modified gravity and cosmological parameters. To install the GitHub version of ISiTGR you can run the following steps in your terminal :

```sh
$ git clone https://github.com/mishakb/ISiTGR
$ cd ISiTGR
Select your compilers by editing the Makefile files inside the principal folder, camb/fortran folder  and source folder.
$ make
```
if there are no errors during compilation, then you should be ready to use ISiTGR.

### Use ISiTGR to calculate cosmological observables
If you want to compute cosmological observables using ISiTGR, then you should install ISiTGR as an extension to CAMB, in order to include MG parameters. 

ISiTGR can be installed using PyPI (https://pypi.org/project/isitgr/) by running
```
$ pip install isitgr [--user]
```
See also https://isitgr.readthedocs.io/en/latest/ for documentation and an example notebook for python-ISiTGR. 

Alternatively, one can use CAMB-ISiTGR by compiling the fortran code:

```sh
$ cd camb/fortran
Select your compilers by editing the Makefile files inside the camb/fortran folder.
$ make
```
if there are no errors during compilation, then you should be ready to use ISiTGR by executing 
```sh
$ ./camb ../inifiles/params_MG.ini
```

In the next table we show some plots to illustrate the different power spectrums that can be produced using the ISiTGR patch for CAMB:

Temperature and lensing potential power spectrum with (μ, η) | TT and TE power Spectrum with (μ, Σ)
:------------------------:|:---------------------:
![](https://dl.boxcloud.com/api/2.0/internal_files/624542586000/versions/662178477600/representations/png_paged_2048x2048/content/1.png?access_token=1!09kbTECN751CPbLqc2tYVzbbKyVf5g2e_Ekk7-ehIGKwduFVwzb8ydWqP1qzcexUoXaVk0ueCU12o7v32Yg_lgqtHu8nvMiIqqbSHAwy5V7Mc9X-It87-xaFCw1zCTg4Zd_JJ0bQUMSCRdJEGbFTqvbRQSqJ_Yqp8bY7jHLqgxyqct1BlPNE1ZjiKiBLmPUyYiOH6k51Y6d_FBqythcERFbMkEM8H8GNmC57XFbQM9FMq_yGmb5_iTj4m54DelnV-P9Eu_22taQoraiNm8WZJkDWh1sxNHwOJmLg7zyQESV1TZsK1w9nX339IZYqxDwKDsC1Y03dnarxDH32XUpy_lWOBBbX5ErxKFONEg2CxgCB121i648mrym5_rMvnIqHaX_S7CYiWFdOrJCweuwjpwGXbMKUPAZOzF788rAOGRAUsdDesqX2MZ2cvhktu8zGUxmdAPyBjVqAPBMkeMlLd5GSGHiuB0B0aYxfqIZyyoMEGcp3TKX26f_zKxCgYSN5XHCSgRVySJiI0Ye76PRAsSqa8njbte_oq9fSYHnU0W2N2we_cgMAjf6hqB4Hk5WboQ..&box_client_name=box-content-preview&box_client_version=2.36.0)   |![](https://dl.boxcloud.com/api/2.0/internal_files/624554338960/versions/662191846960/representations/png_paged_2048x2048/content/1.png?access_token=1!EvBxIyCnPn0hz5vCUFXcnFs3nRV2g0WCprwNcMlYVmChaEGpW0-zo-Sw3BvUO2psuvlmZeJHGbiWklYizw6WoT8X36uDU-d1CO13IX25UWl0pQSdO30Nat0UYmoNBy6_KkIEC9o8H0j0IcSgAcbtlyYM2HgpBoXYaQabFGZMByXE9MyD3woy8pmH8bD9BUjvL8J67o2avddmVfEEyCXPDG46zk_ZHDD4rfmWkpE9tr5MU3RmOqvjT37ZVj4uAx5FULijVuv6qo2d7rVTcq7H8etPOSiTAEgzmyLEQnDCy9F8Im7mwAMtbOm-hbVz-1YaejZm3diyY4F6SPjbUPesn3PonOIIBm-yhjQ6B26aM7y8UVDFUGNGy2mBTUMDcYQUdHmcSkmiikL22a-VK1rSSCCs_hwnY8b0pQM_h74XQHO0MNnaKCJnEGbM69hnfXhjwh2NiXzbddO6EALqC2lJF3bx-QG2S7dvQ1M6sLKZmMazcPmivEUT3N06vAdKy1qfa9-DZtiAKtnmG-sXVSZgxzlD_JENbfV1A4mmFNePy_4HYW23L3RMldmtHS9LTquL_Q..&box_client_name=box-content-preview&box_client_version=2.36.0)  |


## How to run ISiTGR
The ISiTGR patch allows to run CAMB and CosmoMC for MG different models. The next flowchart shows the different parameterizations and methods that ISiTGR is able to work with:
![Parametrizations](https://dl.boxcloud.com/api/2.0/internal_files/624555886276/versions/662193100276/representations/png_paged_2048x2048/content/1.png?access_token=1!hhFT7zg7v3pwu8tz4HHOOoeB3Tgv98otuEqzdtl_AT-_Dye57HLKQINle7F-4ClOw3EKMtiCj93htBLNn0bgl6a47lwxs8mHb9CD1vQjhrkCqz9XzGRdpBg-8x9pnzl-a7dbNE9K29KxFBk8_HPJ4kAgWYlnk13ZoLv8AZGX3I3KyVa9BLp1FDA931mOWZoNCXcVfR7nS4b4g3GQTKZ9war_TygZE44qaV-EU9BAcbhKyWm8i_tiX2t9BlQN-RKd3sSq6f5U7d83uIBLkKUTWId7psy6NrvJ9M1mlVny_frs73W6jVkB1k3oew3GSjiNwCCiq-ZRRUEV3YWCKGCFwpHLF-WIirLwLWV9Dic7FJo3k92fkNgJwNnrO85EQyQhjXGR6Y39a0owETog1gCVJA7JYesENJ-YIVTvbA-wYehrdKJT-HnbRLgSU8z2azPo1lNWhLjtxG6dSTp9soWuDs13Fnyvs_5k82FLorBqsrXOwjhNAyOFDH_EfR8EawVHppDnYL4VCajzqOWp7WMqW_aZFR67Eug9Z5tNrLQHzASlRrrzq03zAluoSBo24YRdTQ..&box_client_name=box-content-preview&box_client_version=2.36.0)
The user can take advantage of the ISiTGR capabilities after compiling the code by modifying the corresponding .ini files. The next flowchart shows the different .ini files that you need to modify in case you want to use the functional form or the binning methods for either CAMB or CosmoMC (you can find further instructions inside each of this files). 
![Implementation](https://dl.boxcloud.com/api/2.0/internal_files/624549348527/versions/662186876927/representations/png_paged_2048x2048/content/1.png?access_token=1!ombitm0Flt-RC1SQOfNF2dcWmzCgixfS7iRI-2i9rwDfiUTLxvtBYSwbJcl8Fyo3DneG5EE7654k6PxgVeQR3_LPGzLWO0WIHzuXIfkSvOV4RlszrrZ2lpxWjqUb97Lk_AB7V6Lpya3jOb0QiJ1-PCzkeyCGkkZAqe2V5s0ODI-mKNKFlRJgDh5-DEZMJsi1eHabwE94JiufATiN08eqXgsbiKRiCkHtg5hytSsBX17Rs6m1IaH2WQ-ILXswHNWr-UcpbmuVsbdV8OOdZMxY_vli2bKsCKc5T74k6MGCERw5hzq9c_B_oRcZPUEItBRB9IXHgbQlrEV3u-fWH-DQB4-o1SBoqBNQgCzsg6dNS-rpDL568RXnr8lygfHMWUioEjupIBWYgo0ydCr_RvuqdXtc7zdMn60KSpRBIg8V8KhiZsGQe6AVj78-a8IjSLgclloLZT2pfMPMaL4fuL2JlRX9MlUfI35xtkNbCwqLM1Dh3I3cuyfKWiwaEcwrahZmPk4PP7UF01_QF55W7CL2dkNzAsnkDoxbr0eWL-siVo5Uea3XvdleRstCBq-zD6J0xA..&box_client_name=box-content-preview&box_client_version=2.36.0)
It is important for the user to remember that the current version of ISiTGR is aimed to work with flat and non-flat models. Moreover, ISiTGR not only implements the contributions of massive neutrinos consistently, but also works with different equations of state for dark energy. The above mentioned is implemented for both functional form and binning methods. If the user wants to run the functional form of ISiTGR, the user needs to run 
```
./cosmomc test_ISiTGR.ini
```
. Otherwise, if the user wants to use the binning methods, then one should run 
```
./cosmomc test_ISiTGR_BIN
``` 

In the next table we show some plots to illustrate some features of the ISiTGR patch for CosmoMC:

Binning Methods |  MG + Curvature  |  Planck 2018 data 
:-------------------------:|:------------------------:|:---------------------:
![](https://dl.boxcloud.com/api/2.0/internal_files/624538642134/versions/662174546934/representations/png_paged_2048x2048/content/1.png?access_token=1!WY0yRBx7TFpouFvdi5KQXKPmfx4SpOwuIkwl023xX4QxXgHKcW_pTLm7xAMliFgB9pwXv9_Q-hUIurQkwsrS8qll8dNjo8ED60cReevJa7Xrqghk82mOfu74huIAcEmQE59PwBaE0PAzlSFdi33b74IhFmswtmfpHmeDSD0tnO7QXwHaR-W4zjhosPozwkUoQ6y5GSdVefdBMk5EC93OeDAtqbGaKjDC2CKfwpJLisT7xX4DLR7duNGo12lUqrOB0W-5q7oQ4SSOIfTkZZeqtj8w8GRElvCkxiKcCEUcC-bgMYk2OghW6X2GkbJdx-FmZKXpY5AdTfm3bKm7PxNmq0p8NIt8v-BaJBKNPenP9rZLJESfR5mvB-mRXUgOMz_ADR3Dv00pkduk1Vrn8t8_5Kl9axEUNZ0_TDvaNDCxZfABJrOOvM2w5FCSKrzr72EjlcmNyQ7foHGyFEkh-P2CcDD5eHoic4zVbfhEFXPo3KoT_ShJTdB8kkkO07v3I5KFcpX7fJPaiGEDfadJq2rjjXCXs2SU2-ZQJ6jNi6K7YM1jcZFe_crqJEQIqTGQ-PBLKA..&box_client_name=box-content-preview&box_client_version=2.36.0)  |  ![](https://dl.boxcloud.com/api/2.0/internal_files/624552947927/versions/662190065927/representations/png_paged_2048x2048/content/1.png?access_token=1!tqW_jyy8N5TGxKtQ3XNXXJz2C1b1z1kSrr0Mf6BgYx7J6JDwWLjixUyW7Dw0qDiGMAA5Dbd38w4EEnRIUMOAsP5fAyVgs9QNJ_JYfb2q7EvDC_Zb5jqGn9oZZDwnj9OrqhO4xv1Fvqiqf5yYeYpK1ddTlNY82S6QRgf8xDhoocKR8XaGEsk7FmpCfsrp81xFdGk89gVY7jeQCVwOzN2iWZg0qNRKntWkN-v81DG1tsyFd6tz-WNwNv0snr8FYfGshb5Lc7S7BSHyvfOP-oTbSSeQUOhJdR1WUElEHpSmYNGnfP22FX9XOMeVHbHHogqWBqDzD5p7h6-smaOOS0WbBRT249I0DaXn2GgI_pGML8yP-hwalcMeiIAb8VCsCpiG7NTAYCumptm0YfiM0aoQe-Iu2dAGJF0ae_VoNwYQQKLCOuZOzFvptSDP3tzticSQJt1sTUyltl8VEXlB19N78vP-y3ABybS9XN-zULDRxxUquxvqNI48JGCMow-sKfGKx54bhchJviKKuDhbhvgJUl3kOpeL2mEY6y5LCClIPfA-hQfECox3IbdtPW84YvAdPQ..&box_client_name=box-content-preview&box_client_version=2.36.0)  |![](https://dl.boxcloud.com/api/2.0/internal_files/624557633619/versions/662195140419/representations/png_paged_2048x2048/content/1.png?access_token=1!aZGl_RXjrnc7w-Wl_GjiZaowNBlhf9pWkBIf13LVHpsZNw6eEv-ZvBlsqg8EiM0sPilc7CCMqg6YYQVo7cAOSYD61Jl5AF7agBp3mTRzPNxmaF7CZnErlOFL8aPZxC7DEbmlMl4-Dh2SvxV1RpBofTcU_8moR4gTcIyVvajTfqKi18bCKnib5lmxV8FLIkvKL79i9Lz2AqkAxJkj_m4Mdv1B-UPV7cCpvX9lfYjAwXaSUh-un5Nr5Oq7RrZqEWi9UFWC3gdTxZ63bgH0pYZbxzJvMSa_vO8gVH0I_cCmtXZrMxMZZftnGwcjHT_15uh0fbje_3D42Ecz3jN3mYXnLx0-p4GhlMfDlMpXLVIJRIq6otiJLxdiqzwZkOJGD4-dPYp7QIQQj2Ap6_MoSGVqnASkI7xKW4gDoLplkHWS4KzxMkWRGsUZNviY3sE6POZXpBdA4PgpuYBg_JZOlTGcRuWBhW1mKRSs2Ntn3EhQKoWxN5R0IUW1F_MSPleOFB2nvxCAWR8XUd6vGfrPm6S-wZ-Vnak1H5bi84SlBV4gpBil3tL6kHzemAIrDkNz-e9QVg..&box_client_name=box-content-preview&box_client_version=2.36.0)  |


## Referencing ISiTGR

We would ask that when using ISiTGR or a modified version of it, you please cite our papers: 
* **ISiTGR: Testing deviations from GR at cosmological scales including dynamical dark energy, massive neutrinos, functional or binned parametrizations, and spatial curvature**: https://arxiv.org/abs/1908.00290  (current ISiTGR version)
by: Cristhian Garcia-Quintero, Mustapha Ishak, Logan Fox and Jason Dossett.
* **Spatial curvature and cosmological tests of general relativity**: https://arxiv.org/abs/1205.2422
by: Jason Dossett and Mustapha Ishak 
* **Testing General Relativity at Cosmological Scales: Implementation and Parameter Correlations**: https://arxiv.org/abs/1109.4583
by: Jason Dossett, Mustapha Ishak and Jacob Moldenhauer
 
as well as the original CAMB paper; the original CosmoMC paper; Additionally please cite the use of any other datasets already included in the original version of CosmoMC.
## Contact

If you have comments, questions, or feedback, please feel free to contact to the contributors of this repository: mishak@utdallas.edu or gqcristhian@utdallas.edu.

<details>
  <summary>Information from previous ISiTGR versions</summary>

-------------------------------------------------------------
# ISiTGR Version 3.01 released in September 2019
We introduce a new version of **I**ntegrated **S**oftware **i**n **T**esting **G**eneral **R**elativity (ISiTGR) which is a patch to the software CAMB and CosmoMC. ISiTGR is intended to test deviations from GR at cosmological scales using available cosmological data sets. While doing so, it allows for various extensions to the standard flat LCDM model. In this new release, we have a combined support for the following:  

* Dynamical dark energy parametrizations with a constant or time-dependent equation of state; 

* A consistent implementation of anisotropic shear to model massive neutrinos throughout the full formalism;

* Multiple commonly used parameterizations of modified growth (MG) parameters; 

* Functional, binned and hybrid time- and scale-dependencies for all MG parameters (expanded from previous version); 

* Spatially flat or curved backgrounds (present in previous version as well). 

The description of the formalism and its implementation in the CMB code, the Integrated Sachs-Wolfe (ISW) effect, and the 3x2 point statistics as well as examples of application to current data sets, can be found in the latest paper on the arXiv. A more technical description of the implementation can be found in the documentation provided in this repository. 

## Installation
The corresponding version of ISiTGR was built based on the CosmoMC July 2018 version. To install the GitHub version of ISiTGR you can run the following steps in your terminal :

```sh
$ git clone https://github.com/mishakb/ISiTGR
$ cd ISiTGR
Select your compilers by editing the Makefile files inside the principal folder, camb folder and source folder.
$ make
$ make isitgr
```
if there are no errors during compilation, then you should be ready to use ISiTGR.

## How to run ISiTGR
The ISiTGR patch allows to run CAMB and CosmoMC for MG different models. The next flowchart shows the different parameterizations and methods that ISiTGR is able to work with:
![Parametrizations](https://drive.google.com/uc?export=view&id=1tgbA9weuGTCAc5qjdnCWkkCgvNu5JVCx)
The user can take advantage of the ISiTGR capabilities after compiling the code by modifying the corresponding .ini files. The next flowchart shows the different .ini files that you need to modify in case you want to use the functional form or the binning methods for either CAMB or CosmoMC (you can find further instructions inside each of this files). 
![Implementation](https://drive.google.com/uc?export=view&id=108szPz_VvQQI9anga74-jcSSUgwmavx1)
It is important for the user to remember that the current version of ISiTGR is aimed to work with flat and non-flat models. Moreover, ISiTGR not only implements the contributions of massive neutrinos consistently, but also works with different equations of state for dark energy. The above mentioned is implemented for both functional form and binning methods. If the user wants to run the functional form of ISiTGR, the user needs to run `./cosmo test_ISiTGR.ini`. Otherwise, if the user wants to use the binning methods, then one should run `./cosmo test_ISiTGR_BIN`. 

In the next table we show some plots to illustrate some features of the ISiTGR patch for CosmoMC:

Binning Methods |  MG + Curvature  |  MG + dark energy + neutrinos 
:-------------------------:|:------------------------:|:---------------------:
![](https://drive.google.com/uc?export=view&id=1deSsgLoSpMYaA6l5Z9-7PTggY1yDgz_O)  |  ![](https://drive.google.com/uc?export=view&id=10BKqA6dDem9rDB1-VnwR1Gcwi0Bu688Y)  |![](https://drive.google.com/uc?export=view&id=1vQxtX3uX73LyL0TgWchsy1uo2FmOkRBY)  |

In the next table we show some plots to illustrate the different power spectrums that can be produced using the ISiTGR patch for CAMB:

Angular power Spectrum for (μ, η) | Angular power Spectrum for (μ, η) allowing scale dependence
:------------------------:|:---------------------:
![](https://drive.google.com/uc?export=view&id=1lB5BRwO5uuTH9EvZvH_C1uUmdrkDdcXw)   |![](https://drive.google.com/uc?export=view&id=1RpZPqidQzZgYK7UWDzLmH_2MlvW0O9ZS)  |


# ISiTGR Version 3.00 released in July 2019
ISiTGR version used in https://arxiv.org/abs/1908.00290 to reproduce Planck 2015 results. The ISiTGR version 3.01 includes minor updates to work with the Planck 2018 data.


# ISiTGR Version 2.01 (previous version information below)

Integrated Software in Testing General Relativity

Version 2.01

Developed by Jason Dossett, Mustapha Ishak, and Jacob Moldenhauer.

What is ISiTGR?

ISiTGR is an integrated set of modified modules for the software package CosmoMC for use in testing whether observational data is consistent with general relativity on cosmological scales. This latest version of the code has been updated to allow for the consideration of non-flat universes. It incorporates modifications to the codes: CAMB, CosmoMC, and the ISW-galaxy cross correlation likelihood code of Ho et al. Also included is our independently developed generalized weak lensing likelihood module with data sets for the CFHTLenS weak lensing tomography of Heymans et al and CFHTLens 2D weak lensing measurements from Kilbinger et al.

A detailed explanation of the modifications made to these codes allowing one to test general relativity are described in our papers: arXiv:1109.4583, arXiv:1205.2422, and arXiv:1501.03119.

How to get ISiTGR

The two versions of ISiTGR have been consolidated into a single package. The three methods of evolving the parameters used to test general relativity, as described in our paper arXiv:1109.4583, are all contained within the code below. Different evolution methods are chosen by using different .ini files and changing options within those files.

This version of ISiTGR is for the July 2015 version of CosmoMC. The original (flat only) verison of ISiTGR as well as builds for other versions of CosmoMC are available here .

Running ISiTGR. (See further below ReadMe-II)
To run ISiTGR you must first download and install CosmoMC (July 2015 version). You can find download instructions, system requirements, and setup instructions for CosmoMC in the CosmoMC ReadMe. Then, simply copy the contents of the unzipped ISiTGR folder to the CosmoMC folder, edit the Makefile for your compilers, compile using "make isitgr", and you are ready to go. To use the functional evolution of the modified gravity parameters, run the code with an ini file based on test_ISiTGR.ini. To use one of the binning methods run using an ini file based on test_ISiTGR_BIN.ini (Please Note: some changes to these ini files may be may be neccessary in order to run with the Planck 2013 likelihoods). Detailed explanations of the specific options available in ISiTGR are availalble in the ISiTGR ReadMe.

Referencing ISiTGR.
We would ask that when using ISiTGR or a modified version of it, you cite: our papers (arXiv:1109.4583 and arXiv:1205.2422), this website; the original CAMB paper; the original CosmoMC paper; and the papers on original ISW-galaxy cross correlation likelihood by Ho et al and Hirata et al. When using any of the CFHTLenS weak lensing data sets, please cite our paper, the corresponding CFHTLens dataset paper (Heymans et al. or Kilbinger et al.), and follow the CFHTLenS publication guidlines given here.  Additionally, please cite the use of any other datasets already included in the original version of CosmoMC.

Version History
2.01 (Released 07/23/15, this version) Updated to July 2015 version of CosmoMC.
2.0 (Released 05/07/15), Major release upgrade: Consolidated the two versions of the code into a single package, major updates to the likelihood modules for compatibility with the Feb. 2015 version of CosmoMC.
1.2 (Released 01/25/15): Updates for Dec. 2013 version of CosmoMC, new CFHTLenS likelihood module.
1.1.1 (Released 05/29/13): Includes bug fixes (Special thanks to Ana Caramete and Lucia Popa).
1.1: Updated code to allow for non-flat universes when varying the MG parameters.
1.02.1 (Released 04/27/12): Fixed bug in CMB_Cls_simple.f90.
1.02: Updated to January 2012 version of CosmoMC.
1.01.1 (Released 04/27/12): Fixed bug in CMB_Cls_simple.f90.
1.01: Updated to October 2011 version of CosmoMC. BAO module now comes included in CosmoMC (Thanks to Antony Lewis).
1.0: Initial release, for August 2011 version of CosmoMC.
Contact
If you have any question or comments about ISiTGR or would like to be updated about with future changes to the code, please feel free to email Mustapha Ishak at: mishak@utdallas.edu and Jason Dossett at: mishak@utdallas.edu.

Readme-II : further info

Introduction

ISiTGR is an integrated set of modified modules for the software package CosmoMC for use in testing whether observational data is consistent with general relativity on cosmological scales. This latest version of the code has been updated to allow for the consideration of non-flat universes. It incorporates modifications to the codes: CAMB, CosmoMC, and the ISW-galaxy cross correlation likelihood code of Ho et al. Also included is our independently developed generalized weak lensing likelihood module with data sets for the CFHTLenS weak lensing tomography of Heymans et al and CFHTLens 2D weak lensing measurements from Kilbinger et al.

To use ISiTGR you must first download and install CosmoMC (July 2015 version). You can find download instructions, system requirements, and setup instructions for CosmoMC in the CosmoMC ReadMe. Then, simply copy the contents of the unzipped ISiTGR folder to the CosmoMC folder, edit the Makefile for your compilers. Compile using the command make isitgr, and you are ready to go.

Basic usage of ISiTGR is the same as CosmoMC. See below for extra options in the .ini files.

Differences between ISiTGR and CosmoMC
Other than the changes to the code to incorporate the parameters used to test general relativy as described in our papers: arXiv:1109.4583 and arXiv:1205.2422, the following additional changes have been made

ini files
There are two input ini files provided with ISiTGR. test_ISiTGR.ini allows you to run the code with the MG parameters evolved using our the functional form evolution. test_ISiTGR_BIN.ini allows you to run the code with the MG parameters evolved using one of the binning methods. The following options are available in boths of these ini files.
To run the ISW-galaxy cross correlation likelihood code, uncomment the line DEFAULT(batch2/ISWHo.ini) in the test_*.ini file that you are using.

To use the CFHTLenS weak lensing tomography data set simply make sure the line DEFAULT(batch2/CFHTLens.ini) is uncommented in test.ini. This file enables use of the full CFHTLenS data sample as described in Heymans et al. If you would like to use one of the other galaxy samples instead you can change that line to:
"Early Type," Red sample - DEFAULT(batch2/CFHTLens_red.ini)
"Late Type," Blue sample - DEFAULT(batch2/CFHTLens_blu.ini)
"Optomized Early Type" sample - DEFAULT(batch2/CFHTLens_rfbb.ini)
These files each contain various settings to enable the different data sets as well as the default parameter ranges for the CFHTLenS intrinsic alignment (IA) nuissance parameter, ACFHTLenS, which allows this parameter to vary from -10 to 10. If you would like to fix this parameter to 0 and thus ignore the (IA) signal you can add the following line at the end of test.ini:
param[CFHTLensA] = 0

To use the CFHTLens 2D weak lensing data set from Kilbinger et al, simply add the line DEFAULT(batch2/CFHTLens_2D.ini) to your test_*.ini.
Additional parameters and options for testing general relativity: Functional Evolution
The file batch2/common_ISiTGR.ini contains lines the line: parameterization = ISiTGR which tells CosmoMC that we want to use the functional evolution of the MG parameters.

The following options have been included in batch1/params_ISiTGR_defaults.ini

The option: Use_R_Function = lets you decide which parameters for testing GR you want to evolve using the funcitonal form shown in our paper. Setting this option to T evolves the parameters Q and R, while setting this option to F evolves Q and Σ (previously called D).
For varying the parameters used to test GR we add the following parameter lines, please refer to our paper for a detailed explanation of each of these parameters.
param[Q0] = 1 0 10 0.05 0.05

param[Qinf] 1

param[Sigma0] = 1 0 10 0.05 0.05

param[Sigmainf] = 1

param[kc] = 0.01

param[s] = 0 0 3 -1 1

The option: Scale_Dependent = lets you decide whether the evolution of the parmaters is scale dependent. Setting this option to T enforces scale dependence, and the parameters ending in inf above should be varied or changed to whatever value you want the parameters to take on small scales.

Additional parameters and options for testing general relativity: Binning Methods
The file batch2/common_ISiTGR_BIN.ini contains lines the line: parameterization = ISiTGR_BIN which tells CosmoMC that we want to use the binning methods to evolve the MG parameters.

The following options have been included in batch2/params_ISiTGR_BIN_defaults.ini

For varying the parameters used to test GR we add the following parameter lines, please refer to our paper for a detailed explanation of each of these parameters.
param[Q1] = 1 0 10 0.05 0.05

param[Q2] = 1 0 10 0.05 0.05

param[Q3] = 1 0 10 0.05 0.05

param[Q4] = 1 0 10 0.05 0.05

param[Sigma1] = 1 0 10 0.05 0.05

param[Sigma2] = 1 0 10 0.05 0.05

param[Sigma3] = 1 0 10 0.05 0.05

param[Sigma4] = 1 0 10 0.05 0.05

param[kc] = 0.01

The option: z_grid_spacing = defines the size of your bins in redshift bins.
The option: Do_Exponential_Binning = lets you choose how you want to transition between scale bins. Setting this option to T sets transitions between scale bins to behave as an exponential function with a decay constant kc (This is the hybrid method described in our paper). Setting this option to F sets transitions between scale bins to behave as a hyberbolic tangent functions (near step functions) with the bins divided at k = kc.

ISW-galaxy cross correlation likelihood code
To implement the ISW-galaxy cross correlations likelihood code of Ho et al, the following files have been added to source/: lrg_2dCl.f90 and iswdata.f90 . The files jl.f90 and lrg_pk.f90 have been removed in this latest version of the code.
The folder bdndz_code which contains c-codes called by the routine Iswlnlike in iswdata.f90 has been added to the root directory. These c-codes are compiled by running the make isitgr command from the top level cosmomc directory.
When Iswlnlike calls the c-program dndz.x it passes input files and output file names to the program. By default, these files are saved in the IO directory within the bdndz_code directory, however you can enable writing of these files to a different directory by setting the environmental variable TMPDIR.

Generalized Weak Lensing Likelihood module and Codes
ISiTGR includes a generalized weak lensing likelihood module, source/WeakLen_Common.f90. This module defines the derived type TCosmologyWLLikelihood and various associated procedures. From this derived type, one can define specific routines for various weak lensing likelihoods. The specific routines relating to the CFHTLens data sets are contained in CFHTLens.f90. All of these files are located in the source directory.
These files are loosely based upon the files contained in the original COSMOS weak lensing likelihood code of Lesgourges et al.

Referencing ISiTGR
We would ask that when using ISiTGR or a modified version of it, you cite: our papers (arXiv:1109.4583 and arXiv:1205.2422), this website; the original CAMB paper; the original CosmoMC paper; and the papers on original ISW-galaxy cross correlation likelihood by Ho et al and Hirata et al. When using any of the CFHTLenS weak lensing data sets, please cite our paper, the corresponding CFHTLens dataset paper (Heymans et al. or Kilbinger et al.), and follow the CFHTLenS publication guidlines given here.  Additionally please cite the use of any other datasets already included in the original version of CosmoMC.

</details> 
