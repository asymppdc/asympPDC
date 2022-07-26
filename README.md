# MATLAB and Octave AsympPDC Toolbox

July 22, 2022

The AsympPDC package is MATLAB/Octave routines and functions comprising to analyze time series or biological signals to determine directional interactions between structures through the Partial Directed Coherence (PDC), which is based on the concept of Granger causality, and the Directed Transfer Function (DTF) in the frequency domain both in three metrics --- Euclidean, diagonal and information --- and rigorous asymptotic statistics providing p-values and confidence interval in the frequency domain. 

## Installation and usage

The AsympPDC package contains MATLAB/Octave compatible mfiles and subfolders you may copy into your local preferred working directory to execute them. To begin with you should run the startup.m script in the MATLAB/Octave command line window to set path and check for the requirements.

```matlab
>> startup
```

In addition to adding the paths, `startup.m` will also check for the presence of the required MATLAB toolboxes (Control, Signal Processing, and Statistical Toolboxes)/Octave packages (Control, Signal, and Statistics). All other routines, including examples and contributed FileExchange files. This is a standalone version. Most likely it will work in the recent versions of Octave --- 6.3.0, 6.4.0 and 7.1.0  (Please report any issue related to compatibility with Octave).

To run all main examples provide in ./examples subdiretory and verify if your installation is most likely working properly, execute

```matlab
>> run_all_examples
```

after `startup.m` script execution. If `run_all_examples.m` completes successfully, it should generate overlapped figures that one can examine, in MATLAB, issueing 

```matlab
>> tilefigs1
```

that will spread the figures on the screen. The tilefigs1 or tilefigs2 function will not work in Octave environment.

## References

The AsympPDC toolbox implementation is based mainly on the following articles and books.

 [1] L.A. Baccala and K. Sameshima. Partial directed coherence: a new concept
     in neural structure determination. Biol Cybern 84:463--474,2001.
     <https://doi.org/10.1007/PL00007990>

 [2] D.Y. Takahashi, L.A.B. Baccala and K. Sameshima, Connectivity inference
     between neural structures via partial directed coherence. J Appl Stat
     34:1259--1273, 2007. <https://doi.org/10.1080/02664760701593065>

 [3] L.A. Baccala, C.S.N. De Brito, D.Y. Takahashi and K. Sameshima. Unified
     asymptotic theory for all partial directed coherence forms. Philos T Roy
     Soc A 371(1997):1--13, 2013. <https://doi.org/10.1098/rsta.2012.0158>
     
 [4] M.J. Kaminski and K.J. Blinowska. A new method of the description of the
    information flow in the brain structures. Biol Cybern 65:203--210,1991.
    <https://doi.org/10.1007/bf00198091>

[5] L.A. Baccala, D.Y. Takahashi and K. Sameshima. Directed transfer
    function: unified asymptotic theory and some of its implications. IEEE T
    Bio-Med Eng 63:2450--2460, 2016. 
    <https://doi.org/10.1109/TBME.2016.2550199>
    
[6] H. Lutkepohl (2005). New Introduction to Multiple Time Series Analysis. 
                         Springer-Verlag, New York. 

[7] S.L. Marple Jr (1987). Digital Spectral Analysis with Application.
                         Prentice-Hall, Englewood-Cliffs. 
                         
[8] T. Schneider and A. Neumaier (2001): Algorithm 808: ARfit - A Matlab package
                         for the estimation of parameters and eigenmodes of
                         multivariate autoregressive models. ACM Trans. Math.
                         Softw., 27:58-–65.

[9] K. Sameshima and L.A. Baccalá Eds. (2014). Methods in Brain Connectivity 
    Inference through Multivariate Time Series Analysis. CRC Press, Boca Raton


## License

These routines are distributed under GNU General Public License v3.0 under
authorship of Koichi Sameshima and Luiz A. Baccalá - July 2022.
