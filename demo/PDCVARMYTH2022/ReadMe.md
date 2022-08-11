# ReadMe

February 1, 2022

The PDCVARMYTH2022 package is MATLAB routines and functions comprising Supplemental Material of the following the article to generate its figures: 

* [Luiz A. Baccalá and Koichi Sameshima (2022) ‘Partial Directed Coherence and the Vector Autoregressive Modelling Myth and
  a Caveat. ](https://www.frontiersin.org/articles/10.3389/fnetp.2022.845327) 

## Installation and usage

The PDCVARMYTH2022 package contains MATLAB mfiles and subfolders you may copy into your local preferred working directory to execute them. To begin with you should run the startup.m script in the MATLAB command line window to set path and check for the requirements.

`> startup`

Other than adding the paths, `startup.m` will also check for the existence of the required MATLAB toolboxes. All other routines, including data simulation and total PDC estimation are included. This is a standalone version. Most likely it will work in the recent versions of Octave 6.3.0, 6.4.0 and 7.1.0 (Please report any issue related to compatibility with Octave).

To run all Examples, choose [1], and [0] otherwise, when requested on `startup.m` script execution.

The figures it generates are similar to those in the [paper](https://www.frontiersin.org/articles/10.3389/fnetp.2022.845327). Differences are due to possibly different random seeds.

To play with the scripts,  you may also run each Example script individually

`> Example1 | Example2 | Example3 | Example4`

## License

These routines are distributed under GNU General Public License v3.0 under
authorship of Koichi Sameshima and Luiz A. Baccalá - January 2022.
