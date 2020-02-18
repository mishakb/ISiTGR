===================
ISiTGR
===================
:ISiTGR: Integrated Software in Testing General Relativity
:Author: Cristhian Garcia Quintero, Mustapha Ishak and Jason Dossett
:Homepage: https://github.com/mishakb/ISiTGR
:Publication: https://arxiv.org/abs/1908.00290  

Description
============

Integrated Software in Testing General Relativity (ISiTGR) which is a patch to the software CAMB and CosmoMC. 
ISiTGR is intended to test deviations from GR at cosmological scales using available cosmological data sets. 
While doing so, it allows for various extensions to the standard flat LCDM model. 
In this new release, we have a combined support for the following:

    Dynamical dark energy parametrizations with a constant or time-dependent equation of state;

    A consistent implementation of anisotropic shear to model massive neutrinos throughout the full formalism;

    Multiple commonly used parameterizations of modified growth (MG) parameters;

    Functional, binned and hybrid time- and scale-dependencies for all MG parameters (expanded from previous version);

    Spatially flat or curved backgrounds (present in previous version as well).

The description of the formalism and its implementation in the CMB code, the Integrated Sachs-Wolfe (ISW) effect,
and the 3x2 point statistics as well as examples of application to current data sets,
can be found in the latest paper on the arXiv. A more technical description of the implementation can be found
in the documentation provided in the ISiTGR repository or the published paper.

See the `ISiTGR python example notebook <https://isitgr.readthedocs.org/en/latest/ISiTGRdemo.html>`_ 
python notebook for illustrative examples of ISiTGR features.

To install the ISiTGR python package on your computer run:

    pip install isitgr [--user]

The --user is optional and only required if you don't have write permission to your main python installation.
You will need the usual dependencies as needed for the CAMB software.

The python wrapper provides a module called "isitgr" to be documented in `Python ISiTGR-CAMB documentation <https://isitgr.readthedocs.io/en/latest/>`_.
