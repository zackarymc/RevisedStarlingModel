% optimzeRRIdriver.m Script that drives revStarOptSelectedRRIfunc. m
%     to optimize selected parameters and plot and output results

% Function revStarOptSelectedRRIfunc is called with inputs and outputs:

% function [Foutputs]=revStarOptSelectedRRI[Finputs] where Foutputs and
%    Finputs are structure comprised as follows:

% Foutputs: Foutputs.tfit,Foutputs.fittedHct,Foutputs.fittedSVs,...
%   Foutputs.fittedParamValues,Foutputs.fittedFuxes,Foutputs.TimeDate

% Finputs:
% Finputs.t,Finputs.RawHct,Finputs.ParamsToOpt,Finputs.ParamICvals,...
%   Finputs.SVICvals,Finputs.HctICval,Finputs.Juf,Finputs.Jinf,...
%   Finputs.patientID,Finputs.LoBnds,Finputs.HiBnds,Finputs.checkfrac,...
%   Finputs.Ntrials,Finputs.fileSaveFlag,Finputs.datafilename,...
%   Finputs.fluxcheckflag

% function [t,FitHct,FitVar,FitParamVals,TimeDate]=revStarOptSelectedRRI(t,...
    % RawHct,ParamsToOpt,ParamICvals,SVICvals,HctICval,Juf,Jinf,patientID,...
    % LoBnds,HiBnds,checkfrac,Ntrials,fileSaveFlag,datafilename,fluxcheckflag)

% inputs: 
% t=time vector,  
% RawHct=Hct(t) data, 
% ParamsToOpt=labels of parameters from [P0,VT0,ChSl,XJvXJs,Pp,phimax,PosSl]
%    to be optimized,
% ParamICvals= initial values for all parameters (those not optimized will 
%    be fixed), 
% SVICvals = initial values for state variables [Vp, Cp, Vi, Ci],
% HctICval= initial Hct value,
% Juf= the ultrafusion time course, 
% Jinf= the time course of any infusion,
% patientID= the patient ID#,
% LoBnds= low bounds for searching (all parameters),
% HiBnds= high bounds for all parameters, 
% checkfrac= fraction of data at end to use to check fit,
% Ntrials= the number of random trials to use,
% fileSaveFlag= flag to save results (0=no, 1=yes),
% datafilename= name of output file, and
% fluxcheckflag= flag to do sum(fluxes) vs. delta volumes check

% outputs:
% t= time vector, same as input,
% FittedHct= fitted hct(t) using optimized parameter values,
% FittedSVs= fitted state variables [Vp Cp Vi Ci],
% FittedParamValues= optimized parameter values for the fitted parameters,
% FittedPressures= hydrostatic and oncotic pressure for optimal model,
% FittedFluxes= Jv, Js and Jl for optimal model,
% TimeDate= time and date stamp when optimization was done, and
% FittedVar= variance of the fit

% ver 5/16/24

clear variables
close all
clc

checkfrac=0.0;% fraction of data at end to ignore in fitting
fluxcheckflag=1;% set to 1 to compute and output all flux computations
fileSaveFlag=1;% if 1, then crate output file and save results
Ntrials=500 %number of randomized starting points to use
noisefactor=0;% add this to create a realistic hct(t)

if fileSaveFlag==1
    diary('Run documentation');
    diary on
end

specialnote='Test Zack scripts and data'

patientnum='RRIpatient250960';
%get input data--this depends on whatever source is used.
% The following is specific for our simulated data (t,hct,Juf,Jinf)
input_title=strcat('Select input data file for subject',num2str(patientnum));
[datafilename,datapath]=uigetfile('*.mat',input_title);% gets pathway and 
% folder name to location of input data file
cd(datapath);%use data in desired data folder
%load .mat patient data
load(datafilename);% loads in t, hct(t), Juf(t) and Jinf(t)
if hct(1)>1
    hct=0.01*hct;
end

rndhct0 = noisefactor*randn(1,length(hct));
rndhct = hct+rndhct0;%Hct with reasonable noise added
dhct=hct-rndhct;
hctvar=var(dhct);%est variance of rndHct
hct=rndhct;

RawHct=hct;

% lines commented out are default values in code to follow
%ParamsToOpt=[{'XJv'},{'XJs'},{'Pp'}];
ParamsToOpt=[{'XJv'},{'XJs'},{'Pp'}];
% ParamICvals=[2,10400,5e-5,1.06,0.2428,20,   38,0.66];%P0,VT0,ChSl,XJv,XJs,Pp,phimax,PosSl
ParamICvals = [2,10400, 5e-5, .15 ,0.41,18.36,38,0.66];%P0,VT0,ChSl,XJv,XJs,Pp,phimax,PosSl
SVICvals=[4000,73.494,17200,24.5955];%VP0,CP0,VI0,CI0
%SVICvals=[4000,60,17200,24.5955];%VP0,CP0,VI0,CI0
HctICval=hct(1);

patientID=patientnum;
LoBnds=[1.9,8000,1e-5,0.4,0.2,15,35,0.56];%P0,VT0,ChSl,XJv,XJs,Pp,phimax,PosSl
HiBnds=[3,12000, 1e-4, 2,  .5,25,55,0.78];
% LoBnds=[1.9,8000,1e-5,0.15,0.15,15,33,0.56];%P0,VT0,ChSl,XJv,XJs,Pp,phimax,PosSl
% HiBnds=[3,12000,1e-4,2,0.4,25,43,0.78];

% Finputs: structure
Finputs.t=t;Finputs.RawHct=RawHct;Finputs.ParamsToOpt=ParamsToOpt;...
Finputs.ParamICvals=ParamICvals;Finputs.SVICvals=SVICvals;...
Finputs.HctICval=HctICval;Finputs.Juf=Juf;Finputs.Jinf=Jinf;...
Finputs.patientID=patientID;Finputs.LoBnds=LoBnds;Finputs.HiBnds=HiBnds;...
Finputs.checkfrac=checkfrac;Finputs.Ntrials=Ntrials;...
Finputs.fileSaveFlag=fileSaveFlag;Finputs.datafilename=datafilename;...
Finputs.fluxcheckflag=fluxcheckflag;

[Foutputs]=revStarOptSelectedRRIfunc(Finputs);% calls optmization process

% Foutputs: structure
tfit=Foutputs.t;fittedHct=Foutputs.fittedHct;
fittedSVs=Foutputs.fittedSVs;fittedParamValues=Foutputs.fittedParamValues;
fittedPressures=Foutputs.fittedPressures;
fittedFluxes=Foutputs.fittedFluxes;TimeDate=Foutputs.TimeDate;
fittedVar=Foutputs.MinVar;