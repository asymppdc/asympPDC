%% asympPDC Package startup script
%
% Run this file to initialize the set up of asympPDC Package.
%
% You may have to customize this script for your computing environment.
%

% (C) Koichi Sameshima & Luiz A. Baccalá, 2022. 
% See file license.txt in installation directory for licensing terms.


%% Set package version

global asymppdc_version;
asymppdc_version.major = 3;
asymppdc_version.minor = 0;

clc
fprintf('                   Initializing asympPDC Package version %d.%d\n', ...
   asymppdc_version.major, asymppdc_version.minor);
fprintf('                   ==========================================\n');
fprintf(1,'\n')
%% Set path

% Add asympPDC root directory and appropriate subdirectories to matlabpath

global asymppdc_root;
%asymppdc_root = fileparts(mfilename('fullpath')); % directory containing this file
asymppdc_root = [pwd '/'];
% essentials

addpath(asymppdc_root);
addpath(fullfile(asymppdc_root,'core'));
addpath(fullfile(asymppdc_root,'extras'));
addpath(fullfile(asymppdc_root,'utils'));
addpath(fullfile(asymppdc_root,'utils','arfit'));
addpath(fullfile(asymppdc_root,'utils','legacy'));
addpath(fullfile(asymppdc_root,'examples'));
addpath(fullfile(asymppdc_root,'examples','models'));

fprintf('[asympPDC startup] Added asymppdc root directory, \n')
fprintf('                      %s, \n',asymppdc_root);
if isOctave()
   warning('off', 'all');
   fprintf('                   and its subdirectories to the OCTAVE search path.\n')
else
   fprintf('                   and its subdirectories to the MATLAB search path.\n')
end
fprintf(1,'\n')
%% Check if Octave
%
if isOctave()
   fprintf('[asympPDC startup] As this is an Octave environment, checking whether control,\n')
   fprintf('                   statistics and signal packages are loaded or not.\n')
   fprintf(1,'\n')
   installed_packages = pkg('list');
   [dummy, N] = size(installed_packages);
   if N ~= 0
      kPkgs = 0;
      pkglist = {'signal', 'statistics', 'control'};  % List of necessary Octave
      % packages.
      for k = 1:N
         for klist = 1:3
            pkgstr = pkglist{klist};
            if strcmpi(installed_packages{k}.name,pkgstr)
               fprintf('                   Package %s is installed and loaded.\n',pkgstr)
               kPkgs = kPkgs + 1;
            end
         end
      end
      fprintf(1,'\n')

      if kPkgs < 3
         fprintf(2,'[asympPDC startup] Warning: some of required packages were not installed or loaded.\n')
         fprintf(2,'                   The required packages are statistics, signal and control. \n')
         fprintf(2,'                   If the package is already installed, load them into Octave \n')
         fprintf(2,'                   environment before proceeding by executing: \n')
         fprintf(2,'                          >> pkg load <package-name>\n')
         fprintf(1,'\n')
         fprintf(2,'                   If not installed, run \n')
         fprintf(2,'                          >> pkg install -forge io\n')
         fprintf(2,'                          >> pkg install -forge <package-name>\n')
         fprintf(1,'\n')
         fprintf(2,'                   See Octave documentation for loading packages on start up \n')
         fprintf(2,'                   editing ".octaverc" file. \n')
      else
         fprintf('[asympPDC startup] All necessary packages are already loaded into Octave environment!\n')
         fprintf('                   You are good to go.\n')

      end
   else
      fprintf(2,'[asympPDC startup] Warning: the control, statistics and signal packages are not load in Octave environment.\n')
      fprintf(2,'                   If packages are already installed, load them into Octave environment before proceeding with analysis\n')
      fprintf(2,'                   through: \n')
      fprintf(2,'                          >> pkg load signal\n')
      fprintf(2,'                          >> pkg load control\n')
      fprintf(2,'                          >> pkg load statistics\n')
      fprintf(2,'                   otherwise install them by executing \n')
      fprintf(2,'                          >> pkg install -forge io\n')
      fprintf(2,'                          >> pkg install -forge control\n')
      fprintf(2,'                          >> pkg install -forge signal\n')
      fprintf(2,'                          >> pkg install -forge statistics\n')
      fprintf(2,'                   then execute "pkg load <package-name> commands. \n')
      fprintf(2,'                   See Octave documentation for loading packages on start up editing ".octaverc" file. \n')
   end
   fprintf(2,'                   If you see any error message concerning availability of package,\n')
   fprintf(2,'                   please take proper action to correct it before continuing.\n')
fprintf(1,'\n')
end

%% Check for dependencies on other MATLAB(R) toolboxes and fileexchange routines
%% needed in asympPDC

% Check if Statistics Toolbox is available - see if ch2cdf is present
if ~isOctave()
   kPkgs = 0;
   if ~isempty(which('chi2cdf')) && ~isempty(which('chi2inv')) ...
         && ~isempty(which('norminv'))
      fprintf('[asympPDC startup] Statistics Toolbox(TM) or Statistics and Machine Learning Toolbox(TM) \n');
      fprintf('                   seems to be present.\n');
      kPkgs = kPkgs +1;
   else
      fprintf(2,'[asympPDC startup] WARNING: Statistics Toolbox(TM) or Statistics and Machine \n');
      fprintf(2,'[asympPDC startup]          Learning Toolbox(TM) does not seem to be present.\n');
      fprintf(2,'[asympPDC startup]          Will not be able to perform asymptotic statistical \n');
      fprintf(2,'[asympPDC startup]          tests for PDC and DTF estimates. \n')
   end
   fprintf(1,'\n')

   % Check if you have Signal Processing Toolbox - check whether pwelch is present

   if ~isempty(which('xcorr')) && ~isempty(which('norminv'))
      fprintf('[asympPDC startup] Signal Processing Toolbox(TM) seems to be present.\n');
      kPkgs = kPkgs +1;
   else
      fprintf(2,'[asympPDC startup] WARNING: Signal Processing Toolbox(TM) does not seem to be present.\n');
      fprintf(2,'[asympPDC startup]          Will not be able to perform mvarresidue routines \n');
      fprintf(2,'[asympPDC startup]          for the residues test for whiteness. \n')
   end
   fprintf(1,'\n')
end

if ~isempty(which('dlyap')) && ~isOctave()
   fprintf('[asympPDC startup] Control System Toolbox(TM) seems to be present.\n');
   kPkgs = kPkgs +1;
else
   fprintf(2,'[asympPDC startup] WARNING: Control System Toolbox(TM) does not seem to be present.\n\n');
   fprintf(  '[asympPDC startup] The lyap.m legacy code from 1986 Control System Toolbox(TM) will be used.\n');
   addpath(fullfile(asymppdc_root,'utils','legacy'));
   fprintf(  '[asympPDC startup] Adding ./code/utils/legacy directory to MATLAB path with lyap.m.\n');
end

if ~isempty(which('lyap'))
   % Testing lyap.m code using Example 1 from modern lyap help.
   %The A matrix is stable, and the Q matrix is positive definite.

   % Define an absolute tolerance
   tol = 1e-10;
   % Three input arguments example
   A = 5;
   B = [4 3; 4 3];
   C = [2 1];
   X = lyap(A,B,C);
   % If correct, these commands should return the following X matrix:
   Xcorrect = [   -0.2000   -0.0500];

   if sum(abs(Xcorrect(:) - X(:))) <= tol
      fprintf('                   The *lyap* function is present and working properly.\n');
   else
      error('Problem with lyap.m function. Verify lyap.m function.')
   end
   clear A B C X Q Xcorrect tol
else
   fprintf(2,'[asympPDC startup] WARNING: Control System Toolbox(TM) not installed,\n');
   fprintf(2,'[asympPDC startup]          and lyap function is not available.\n');
end
fprintf(1,'\n')

% Check if the ARFit package is present

if ~isempty(which('arfit')) && ~isempty(which('arfitcaps'))
   fprintf('[asympPDC startup] ARFit Toolbox(TM) seems to be present.\n');
else
   fprintf(2,'[asympPDC startup] WARNING: ARFit(TM) does not seem to be present.\n');
   fprintf(2,'[asympPDC startup]          Will not be able to perform asymptotic statistics \n');
   fprintf(2,'[asympPDC startup]          analyses for both PDC and DTF. \n')
end
fprintf(1,'\n')


if ~isempty(which('suplabel')) && ~isempty(which('suptitle'))
   if isOctave()
      fprintf('[asympPDC startup] The xplot function will use suptitle.m (Drea Thomas, 1995).\n\n');
      fprintf('               (*) As suplabel.m (Ben Barrowes, 2004) does not work in Octave,\n');
      fprintf('                   the y- and x-axis labels will be missing in the matrix layout subplots, .\n');

   else
      fprintf('[asympPDC startup] The xplot function will use suplabel.m (Ben Barrowes, 2004), and\n');
      fprintf('                                               suptitle.m (Drea Thomas, 1995).\n');
   end
   fprintf(1,'\n')
if ~isempty(which('shadedplot')) && ~isempty(which('boundedline'))
   fprintf('                   And, also, if required,  shadedplot.m (Dave Van Tol, 2008), and\n');
   fprintf('                                            boundedline.m (Kelly Kearney, 2012).\n');
   fprintf(1,'\n')
end
end

if ~isOctave()
   if ~isempty(which('tilefigs')) || ~isempty(which('tilefigs1')) || ~isempty(which('tilefigs2'))
      fprintf('[asympPDC startup] You may quickly visualize multiple figure windows issuing either\n');
      fprintf('                   tilefigs1 (Peter J Acklam,2003) or tilefigs2 (Brendan Tracey, 2012) command.\n');
   end
end
fprintf(1,'\n')

%% Initialize default random number stream
% Initialize rng to avoid predictability of sessions

rng('default')
rng('shuffle')
fprintf('[asympPDC startup] Random number generator initialized using MATLAB recommended procedure.\n');
fprintf(1,'\n')

%% Enable all warnings

if isOctave()
   warning('off', 'all');
else
   warning on all
end

if kPkgs < 3
   if isOctave()
      fprintf(1,'\n')
      fprintf(2,'[asympPDC startup] WARNING: One or more required Octave package(s) are missing.\n');
      fprintf(2,'[asympPDC startup]          Install the packages before proceeding. \n');
      fprintf(2,'[asympPDC startup]          (You may re-run ''startup.m'').\n\n\n')
   else
      fprintf(1,'\n')
      fprintf(2,'[asympPDC startup] WARNING: One or more required MATLAB toolboxes are missing.\n');
      fprintf(2,'[asympPDC startup]          Statistics, Control Systen and Signal Processing Toolboxes \n');
      fprintf(2,'[asympPDC startup]          are necessary for running asympPDC Package. \n');
      fprintf(2,'[asympPDC startup]          (You may re-run ''startup.m'').\n\n\n')
      
   end
else
   fprintf('[asympPDC startup] All warnings enabled\n');
   fprintf(1,'\n')
   fprintf('[asympPDC startup] Initialization completed (you may re-run ''startup.m'' if necessary).\n\n\n');
end
% <startup.html back to top>

