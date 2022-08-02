# MATLAB and Octave AsympPDC Toolbox

July 22, 2022

The AsympPDC package is MATLAB/Octave routines and functions for the analysis of time series or biological signals to infer directional interactions between structures through the Partial Directed Coherence (PDC), which is based on the concept of Granger causality, and the Directed Transfer Function (DTF) in the frequency domain both in three metrics --- Euclidean, diagonal and information --- and rigorous asymptotic statistics providing p-values and confidence interval in the frequency domain. 

## Installation and usage

The AsympPDC package contains MATLAB/Octave mfiles and subfolders that you may copy into your local preferred working directory to execute them. To begin with you should go to the package directory and run the startup.m script in the MATLAB/Octave command line window that will set paths and check for the requirements.

```matlab
>> startup
```

In addition to adding the paths, `startup.m` will also check for the presence of the required MATLAB toolboxes (Control, Signal Processing, and Statistical Toolboxes) or Octave packages (Control, Signal, and Statistics). This is a standalone version. Most likely it will work in the recent versions of Octave --- 6.3.0, 6.4.0 and 7.1.0  (Please report or suggest correction to any issue related to compatibility with Octave).

To run all main examples provided in ./examples subdiretory and verify if your installation is working properly, execute

```matlab
>> run_all_examples
```

If `run_all_examples.m` completes successfully, congratulation, you should see 28 overlapped figures that could be examine, in MATLAB, issueing 

```matlab
>> tilefigs1
```

that will spread the figures on the screen. The tilefigs1 or tilefigs2 function does not work so far in the Octave environment.



## Our diagrammatic view of connectivity measures

The measures inside the yellow area have been implemented in the asympPDC packages



![](./connectivity_measures_in_asymppdc.png)

## References

The AsympPDC toolbox implementation is based mainly on the following articles and books:

 [1] L.A. Baccala and K. Sameshima (2001a). Partial directed coherence: a new concept
     in neural structure determination. *Biol Cybern* **84**:463--474.
     <https://doi.org/10.1007/PL00007990>

 [2] D.Y. Takahashi, L.A.B. Baccala and K. Sameshima (2007), Connectivity inference
     between neural structures via partial directed coherence. *J Appl Stat*
     **34**:1259--1273. <https://doi.org/10.1080/02664760701593065>

 [3] L.A. Baccala, C.S.N. De Brito, D.Y. Takahashi and K. Sameshima (2013). Unified
     asymptotic theory for all partial directed coherence forms. *Philos T Roy
     Soc A* **371**:1--13. <https://doi.org/10.1098/rsta.2012.0158>

 [4] M.J. Kaminski and K.J. Blinowska (1991). A new method of the description of the
    information flow in the brain structures. *Biol Cybern* **65**:203--210.
    <https://doi.org/10.1007/bf00198091>

[5] L.A. Baccala, D.Y. Takahashi and K. Sameshima (2016). Directed transfer
    function: unified asymptotic theory and some of its implications. *IEEE T
    Bio-Med Eng* **63**:2450--2460. 
    <https://doi.org/10.1109/TBME.2016.2550199>

[6] H. Lutkepohl (2005). New Introduction to Multiple Time Series Analysis. 
                         Springer-Verlag, New York. 

[7] S.L. Marple Jr (1987). Digital Spectral Analysis with Application.
                         Prentice-Hall, Englewood-Cliffs. 

[8] T. Schneider and A. Neumaier (2001). Algorithm 808: ARfit - A Matlab package
                         for the estimation of parameters and eigenmodes of
                         multivariate autoregressive models. *ACM Trans Math
                         Softw* **27**:58-–65. <https://doi.org/10.1145/382043.382316>

[9] K. Sameshima and L.A. Baccalá Eds. (2014). Methods in Brain Connectivity 
    Inference through Multivariate Time Series Analysis. CRC Press, Boca Raton.
    <https://doi.org/10.1201/b16550>

## License

These routines are distributed under GNU General Public License v3.0 under
authorship of Koichi Sameshima and Luiz A. Baccalá - July 2022.
