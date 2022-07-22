# MATLAB and Octave AsympPDC Toolbox

> July 22, 2022

The AsympPDC package is MATLAB/Octave routines and functions comprising to analyze time series or biological signals to determine directional interactions between structures through Partial Directed Coherence (PDC) , which is based on the concept of Granger causality, and Directed Transfer Function (DTF) in the frequency domain in their three metrics --- Euclidean, diagonal and informational --- and rigorous asymptotic statistics providing p-values and confidence interval also in the frequency domain. The AsympPDC toolbox implementation is based on the following main articles and books.


## Installation and usage

The AsympPDC package contains MATLAB/Octave compatible mfiles and subfolders you may copy into your local preferred working directory to execute them. To begin with you should run the startup.m script in the MATLAB/Octave command line window to set path and check for the requirements.

`> startup`

In addition to adding the paths, `startup.m` will also check for the presence of the required MATLAB toolboxes (Control, Signal Processing, and Statistical Toolboxes)/Octave packages (Control, Signal, and Statistics). All other routines, including examples and contributed FileExchange files. This is a standalone version. Most likely it will work in the recent versions of Octave --- 6.3.0, 6.4.0 and 7.1.0 (Please report any issue related to compatibility with Octave ).

To run all main examples provide in ./examples subdiretory to verify if your installation is working properly, execute

```matlab
>> run_all_examples
```

after `startup.m` script execution. If `run_all_examples.m` completes successfully, it should generate overlapped figures that one can examine, in MATLAB, issueing

```matlab
>> tilefigs1
```

that will spread the figures on the screen. This function does not work in Octave environment.

## License

These routines are distributed under GNU General Public License v3.0 under
authorship of Koichi Sameshima and Luiz A. Baccalá - July 2022.
