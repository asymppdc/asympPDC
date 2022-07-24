%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sunspot and melanoma data 1936 - 1972
%     sunspot --> melanoma 

more off
warning off

tic;

andrews_herzberg
shg
pause(2)

%==========================================================================
% Baccala & Sameshima. Partial directed coherence: a new concept in neural
% structure determination. Biol. Cybern. 84:463-474, 2001.
% 
% Example VAR(5) with loop and feedback 4<=>5 Example 3
baccala2001a_ex3
shg
pause(2)

% Example VAR(5) with loop and feedback  (4<=>5) Example 4
baccala2001a_ex4
shg
pause(2)

% Example VAR(5) with loop and feedback
baccala2001a_ex5
shg
pause(2)

%==========================================================================
% Baccala, L. A. & Sameshima, K. (2001b) Overcoming the limitations of
% correlation analysis for many simultaneously processed neural structures.
% Progress in Brain Research, 130, pp. 33--47. 
%
% Example Model I - (VAR(5) + feedback) + independent VAR(2)

baccala2001b_model1_feedback
shg
pause(2)

% Example Model II Left and Right Hemispheres EEG emulation
baccala2001b_model2
pause(2)

% Example Modified Model II Left and Right Hemispheres EEG emulation

baccala2001b_model2_variant
pause(2)

% Eichler. On the evaluation of information flow in multivariate systems 
% by the directed transfer function.
%            Biol Cybern (2006) 94: 469?482
%  Example - Three-dimensional VAR[2].
eichler2006_ex1              
pause(2)

% Eichler. On the evaluation of information flow in multivariate systems 
% by the directed transfer function.
%            Biol Cybern (2006) 94: 469?482
%
%  Example - two-dimensional VAR[4].
eichler2006_ex2              
pause(2)

%==========================================================================
% Gourevitch, Bouquin-Jeannes, Faucon
% Linear and nonlinear casuality between signals: methods, examples and
% neurophysiological applications. Biol Cybern 95:349-369, 2006.
%==========================================================================

% Example Model 2: Linear bivariate model with bidirectional influence 
% with common source
gourevitch2006_model2
shg
pause(2)

%==========================================================================
% Guo, Wu, Ding & Feng. Uncovering interactions in frequency domains. PLoS
% Computational Biology, 4(5):1-10, February 8, 2008. 
%==========================================================================
% Page 2 Toy Model Example VAR(5) + VAR(2) with loop and feedback
guo2008_linear
shg
pause(2)

%==========================================================================
% Schelter, Winterhalder, Eichler, Peifer,Hellwig, Guschlbauer, Lucking,
% Dahlhaus, Timmer. Testing for directed influences among neural 
% signals using partial directed coherence. Journal of Neuroscience Methods
% 152:210-218, 2005.
%==========================================================================
% Example VAR(4) 
schelter2005
shg
pause(2)

%==========================================================================
% J Physiology - Paris 99:37-46, 2006.
% Direct or indirect? Graphical models for neural oscillators
% Bjorn Schelter, Matthias Winterhalder, Bernhard Hellwig, 
% Brigitte Guschlbauer, Carl Hermann Lucking, Jens Timmer
%==========================================================================
% 
% Example VAR(4) 
schelter2006
shg
pause(2)

%==========================================================================
% Schelter, Timmer, Eichler. Assessing the strength of directed influences
% among neural signals using renormalized PDC. J Neurosci Methods (2009 in
% press). 
%==========================================================================
% 3.1 Vector autoregressive process I (pag. )
schelter2009_vap1
shg
pause(2)
schelter2009_vap2            
pause(5)

%==========================================================================
% Winterhalder et al.Comparison of linear signal processing techniques to
% infer directed interactions in multivariate neural systems. 
% Signal Processing 85:2137--2160, 2005.
%==========================================================================
% Variant of Random Independent Process with 7 variables')
% sigma1=500; sigma2=1; sigma3=500; sigma4=1; sigma5=1; sigma6=1; sigma7=1;
%
% Example: Seven random independent variables 

winterhalder2005_variant

disp(repmat('=',1,100))
disp(['=============== ALL EXAMPLES SUCCESSFULLY FINISHED ' repmat('=',1,49)])
disp(repmat('=',1,100))

% fbaccala2001a_ex3            
% fbaccala2001a_ex4            
% fbaccala2001a_ex5            
% fbaccala2001b                
% fgourevitch2006              
% fguo2008_linear              
% fschelter2005                
% fschelter2006                
% fschelter2009_vap1           
% fschelter2009_vap2           
% fwinterhalder2005_variant    
% gourevitch2006               
% guo2008_linear               
% run_all_examples             
% schelter2005                 
% schelter2006                 
     

elapsedtime = toc