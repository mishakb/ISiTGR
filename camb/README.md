[![N|Solid](https://www.utdallas.edu/~jdossett/images/banner_bkgd_isitgr3.jpg)](https://www.utdallas.edu/~jdossett/images/banner_bkgd_isitgr3.jpg)

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
Select your compilers editing the Makefile files.
$ make
$ make isitgr
```
if there are no errors during compilation, then you should be ready to use ISiTGR.

## How to run ISiTGR
The ISiTGR patch allows to run CAMB and CosmoMC for MG different models. In case you want to produce different kinds of power spectrum using CAMB and ISiTGR, the next flowchart shows the different parameterizations and methods that ISiTGR is able to work with:
[![N|Solid](https://drive.google.com/uc?export=view&id=1_jDEFNdN_K9K-i5Ajaz12a5t6U5sWhEI)](https://drive.google.com/uc?export=view&id=1_jDEFNdN_K9K-i5Ajaz12a5t6U5sWhEI)
It is important for the user to remember that the current version of ISiTGR is aimed to work with flat and non-flat models. Moreover, ISiTGR not only implements the contributions of massive neutrinos in a consistent way, but also works with different equations of state for dark energy. The above mentioned is implemented for both Functional Form and Binning methods. If the user wants to run the binning form of ISiTGR, then the user needs to run modify the file `params.ini` and set Binning_Method = F. Otherwise, the user will be using the functional forms for the MG parameters. In either case, the user has to modify either the file `params_ISiTGR.ini` or the file `params_ISiTGR_BIN.ini`, for functional form and binning methods, respectively. 

In the next table we show some plots to illustrate of the different power spectrum that can be produced using ISiTGR:

Angular power Spectrum for (μ, η) | Angular power Spectrum for (μ, η) allowing scale dependence
:------------------------:|:---------------------:
![](https://drive.google.com/uc?export=view&id=1lB5BRwO5uuTH9EvZvH_C1uUmdrkDdcXw)   |![](https://drive.google.com/uc?export=view&id=1RpZPqidQzZgYK7UWDzLmH_2MlvW0O9ZS)  |

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
# ISiTGR Version 3.00 released in July 2019
ISiTGR version used in https://arxiv.org/abs/1908.00290 to reproduce Planck 2015 results. The ISiTGR version 3.01 includes minor updates to work with the Planck 2018 data.

</details> 