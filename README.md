[![N|Solid](https://www.utdallas.edu/~jdossett/images/banner_bkgd_isitgr3.jpg)](https://www.utdallas.edu/~jdossett/images/banner_bkgd_isitgr3.jpg)

# ISiTGR Version 3.1 released in February 2020 (with python wrapper)
We introduce a new version of **I**ntegrated **S**oftware **i**n **T**esting **G**eneral **R**elativity (ISiTGR) which is a patch to the software CAMB and CosmoMC. ISiTGR is intended to test deviations from GR at cosmological scales using available cosmological data sets. While doing so, it allows for various extensions to the standard flat LCDM model. In this new release, we have a combined support for the following:  

* Dynamical dark energy parametrizations with a constant or time-dependent equation of state; 

* A consistent implementation of anisotropic shear to model massive neutrinos throughout the full formalism;

* Multiple commonly used parameterizations of modified growth (MG) parameters; 

* Functional, binned and hybrid time- and scale-dependencies for all MG parameters (expanded from previous version); 

* Spatially flat or curved backgrounds (present in previous version as well). 

Additionally, **this new version of ISiTGR provides a python wrapper**. The python wrapper extends the original CAMB's python wrapper to work with the different MG parameterizations provided in ISiTGR, allowing the user to obtain the power spectra and transfer functions including phenomenological modified gravity parameters.

The description of the formalism and its implementation in the CMB code, the Integrated Sachs-Wolfe (ISW) effect, and the 3x2 point statistics as well as examples of application to current data sets, can be found in the latest paper on the arXiv. A more technical description of the implementation can be found in the documentation provided in this repository. 

## Installation (to run MCMC)
The corresponding version of ISiTGR was built based on the CosmoMC July 2019 version. If you install ISiTGR from the principal folder, you will be allowed to do MCMC sampling to constraint modified gravity and cosmological parameters. To install ISiTGR from the GitHub repository you can run the following steps in your terminal :

```sh
$ git clone https://github.com/mishakb/ISiTGR
$ cd ISiTGR
Select your compilers by editing the Makefile files inside the principal folder, camb/fortran folder  and source folder.
$ make
```
if there are no errors during compilation, then you should be ready to use ISiTGR.

### Use ISiTGR to calculate cosmological observables
If you want to compute cosmological observables using ISiTGR, then you should install ISiTGR as an extension to CAMB, in order to include MG parameters. 

In this new version, ISiTGR can be installed using PyPI (https://pypi.org/project/isitgr/) by running
```
$ pip install isitgr [--user]
```
See also https://isitgr.readthedocs.io/en/latest/ for documentation and an example notebook for python-ISiTGR. 

Alternatively, one can use ISiTGR by compiling the fortran code:

```sh
$ cd camb/fortran
Select your compilers by editing the Makefile files inside the camb/fortran folder.
$ make
```
if there are no errors during compilation, then you should be ready to use the fortran ISiTGR code by executing 
```sh
$ ./camb ../inifiles/params_MG.ini
```

In the next table we show some plots to illustrate the different power spectrums that can be produced using the ISiTGR patch for CAMB:

Temperature and lensing potential power spectrum with (μ, η) | TT and TE power Spectrum with (μ, Σ)
:------------------------:|:---------------------:
![](https://imagizer.imageshack.com/img922/2492/OSyYNQ.png)   |![](https://imagizer.imageshack.com/img924/6177/MJeYME.png)  |


## How to run ISiTGR
The ISiTGR patch allows to run CAMB and CosmoMC for MG different models. The next flowchart shows the different parameterizations and methods that ISiTGR is able to work with:
![Parametrizations](https://imagizer.imageshack.com/img922/5735/VwisJk.png)
The user can take advantage of the ISiTGR capabilities after compiling the code by modifying the corresponding .ini files. The next flowchart shows the different .ini files that you need to modify in case you want to use the functional form or the binning methods for either CAMB or CosmoMC (you can find further instructions inside each of this files). 
![Implementation](https://imagizer.imageshack.com/img921/710/TDBMMJ.png)
It is important for the user to remember that the current version of ISiTGR is aimed to work with flat and non-flat models. Moreover, ISiTGR not only implements the contributions of massive neutrinos consistently, but also works with different equations of state for dark energy. The above mentioned is implemented for both functional form and binning methods. If the user wants to run the functional form of ISiTGR, the user needs to run 
```
./cosmomc test_ISiTGR.ini
```
Otherwise, if the you want to use binning methods, instead you should run 
```
./cosmomc test_ISiTGR_BIN.ini
``` 

In the next table we show some plots to illustrate some features of the ISiTGR patch for CosmoMC:

Binning Methods |  MG + Curvature  |  Planck 2018 data 
:-------------------------:|:------------------------:|:---------------------:
![](https://imagizer.imageshack.com/img921/4875/UF2JJq.png)  |  ![](https://imagizer.imageshack.com/img922/3460/BnmcTM.png)  |![](https://imagizer.imageshack.com/img922/4962/o1LWlk.png)  |


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
