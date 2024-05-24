% function version of revstarlingoptselectedLHS_RRI.m This uses the revised  
% Starling model based on the de los Reyes 2016 model with the substitution  
% of glycocalyx dynamics of capillary exchange using Facchini, Chapple and 
% Possenti methods to optimize selected parameters. Uses fminsearchbnd
% with revstarlingmodel.m to compute model
%
% 3/13/24: functional version of optimization routine 
% 3/17/24: this version does a final fminsearchbnd using the parameter
% set with the minimum variance from the Ntrial loop as the initial values
% 4/14/24: New version uses structures to pass inputs and results. This
% function is called from optimzeRRIdriver.m

function [Fresults]=revStarOptSelectedRRIfunc(Finputs)

% inputs: Finputs structure
% t=time vector,  
% RawHct=Hct(t) data, 
% ParamsToOpt=labels of parameters from [P0,VT0,ChSl,XJvXJs,Pp,phimax,PosSl]
%    to be optimized,
% ParamICvals= initial values for all parameters (those not optimized will 
%    be fixed), 
% SVICvals = initial values for state variables [Vp, Cp, Vi, Ci],
% HctICval= initial Hct value,
% Juf= the ultrafusion time course (negative since always out of Vp), 
% Jinf= the time course of any infusion (positive since always into Vp),
% patientID= the patient ID#,
% LoBnds= low bounds for searching (all parameters),
% HiBnds= high bounds for all parameters, 
% checkfrac= fraction of data at end to use to check fit,
% Ntrials= the number of random trials to use. 
% fileSaveFlag= flag to save results (0=no, 1=yes), 
% fluxcheckflag= flag to compute and display fluxes and sums, and 
% datafilename= name of input data file.

% outputs: Fresults structure
% t= time vector, same as input,
% FittedHct= fitted hct(t) using optimized parameter values,
% FittedVar= variance of the fit,
% FittedParamValues= optimized parameter values for the fitted parameters,
% FittedPressures= hydrostatic and oncotic pressures for optimal model,
% FittedFluxes= fluxes calculated from optimal parameter values, and
% FittedTimeDate= time and date stamp when optimization was done
% FittedMinVar= minimum variance of optimal model

TimeDate=char(datetime);

t=Finputs.t;RawHct=Finputs.RawHct;ParamsToOpt=Finputs.ParamsToOpt;
ParamICvals=Finputs.ParamICvals;SVICvals=Finputs.SVICvals;
HctICval=Finputs.HctICval;Juf=Finputs.Juf;Jinf=Finputs.Jinf;
patientID=Finputs.patientID;LoBnds=Finputs.LoBnds;HiBnds=Finputs.HiBnds;
checkfrac=Finputs.checkfrac;Ntrials=Finputs.Ntrials;
fileSaveFlag=Finputs.fileSaveFlag;fluxcheckflag=Finputs.fluxcheckflag;
datafilename=Finputs.datafilename

solverinfo='solver: ode15s';

MinMethod='fminsearchbnd';
approachinfo=['Approach: Revised Starling model using ',MinMethod];

allParamNames={'P0', 'VT0', 'ChSl', 'XJv', 'XJs', 'Pp', 'phimax', 'PosSl'};

% pick out parameter ids to optimize
paramidsToOpt=[];
paramidsToFix=[];
j=1;
for i=1:8% identify parameters to be optimized
    if strcmp(ParamsToOpt(j),allParamNames(i))==1
        paramidsToOpt(j)=i;
        if j<length(ParamsToOpt)
            j=j+1;
        end
    else
        paramidsToFix(length(paramidsToFix)+1)=i;
    end
end

thetaNames=allParamNames(paramidsToOpt) 
nparams2opt=length(paramidsToOpt);

ModelID=['revStarOptSelectedRRI: optimizing parameters: ',num2str(paramidsToOpt)];

if fileSaveFlag==1 %set up output file
    % construct unique data file name
    pn=char(datafilename);
    dotlocus=find(pn=='.');
    subjid=string(pn(1:dotlocus-1));
    datestr=char(datetime('now','TimeZone','local','Format','yy-MM-dd_HH-mm-ss'));
    outfilename=strcat(subjid,{'_'},datestr);
    xlfname=strcat(outfilename,'xlsx');%Excel output file
    
    %setup data folder for this patient and data set
    resultfoldername=strcat(outfilename,'results');

    selpath=uigetdir();
    cd(selpath)
    mkdir(char(resultfoldername));
    cd(char(resultfoldername));% this is where results will be stored

    diary 'Run documentation.txt';
    diary on
    ModelID
    Ntrials
    patientID
end

%% Check data and set time limits
% Input data format checking
%  check to make sure hct is a row vector

hct=RawHct;% entire hct data vector
if size(hct,1)~=1 %row vector is (n,1)
    hct=hct';%convert column to row vector
end

if hct(1)>1% check to see if hct is a fraction
    hct=.01*hct; % correct to fraction
end

tall=t;%entire data set including what will be fitted and that saved for checking
jufall=Juf;
hctall=hct;
jinfall=Jinf;

% Ensure that all data arays are row vectors (1xn)
[~,ncols]=size(tall);
if ncols==1 %column vector needs to be transposed
    tall=tall';
end
[~,ncols]=size(hctall);
if ncols==1 %column vector needs to be transposed
    hctall=hctall';
end
[~,ncols]=size(jufall);
if ncols==1 %column vector needs to be transposed
    jufall=jufall';
end

deltat=tall(2)-tall(1);%sampling interval for data
tspanall=[tall(1) tall(end)];%time limits for entire (fitted +checked) system solution
ndatatotal=length(tall);% number of data points total

%checkfrac=fraction of end of data to use for checking
%   do not fit this fraction at end of data
fitfrac=(1-checkfrac);% fraction of data to use for fitting
%select data for fitted region
tlow=0;thigh=t(fix(fitfrac*ndatatotal));%time limits in sec for data span to fit
nlow=find(tall>=tlow,1);%index of lowest data point to use
nhigh=find(tall>=thigh,1);%index of highest data point to use
ndatafit=nhigh-nlow+1;% number of data points to use for fitting

%select data for fitted region
tfit=tall(nlow:nhigh);%times for fitted data
hctfit=hctall(nlow:nhigh);%hct values for fitted data
juffit=jufall(nlow:nhigh);%Juf values for fitted data
tspanfit=[tfit(1) tfit(end)]; %start and end times span fitted data

if ndatafit==ndatatotal% all data used for fit, none for checking
    tcheck=[];
    hctcheck=[];
    jufcheck=[];
    tspancheck=[];
else %values for checked data
    tcheck=tall(nhigh+1:end);
    hctcheck=hctall(nhigh+1:end);
    jufcheck=jufall(nhigh+1:end);
    tspancheck=[tcheck(1) tcheck(end)];% time span for checked data
end

%select data for low and high check regions
if ndatafit==ndatatotal% all data used for fit, none for checking
    tchecklow=[];
    hctchecklow=[];
    jufchecklow=[];
    tspanchecklow=[];
    tcheckhigh=[];
    hctcheckhigh=[];
    jufcheckhigh=[];
    tspancheckhigh=[];
else %data at end not included in fitting-used for checking
    tchecklow=tall(1:nlow);%times below fitted data span
    hctchecklow=hctall(1:nlow);
    jufchecklow=jufall(1:nlow);
    tspanchecklow=[tchecklow(1) tchecklow(end)];% time span for checked data
    nchecklow=length(tchecklow);%number of data points prior to fitted range

    tcheckhigh=tall(nhigh:end);%time span above fitted data span, used to check fit
    hctcheckhigh=hctall(nhigh:end);
    jufcheckhigh=jufall(nhigh:end);
    tspancheckhigh=[tcheckhigh(1) tcheckhigh(end)];% time span for checked data
    ncheckhigh=length(tcheckhigh);%number of data points prior to fitted range
end

%Applied flux = sum of Juf out of and Jinf into plasma
appliedFlux=-Juf+Jinf;%net flux into Vp
% sets up interpolation grid for subsequent use by ode solver
interpApplFlux=griddedInterpolant(t,appliedFlux,'nearest');

%State variable and Hct initial conditions for data
Hct0=HctICval;
Vp0=SVICvals(1);
Cp0=SVICvals(2);
Vi0=SVICvals(3);
Ci0=SVICvals(4);

% set baseline parameter values- default values
% these will be used for those parameters not in the to-be-optimized set
P0init=ParamICvals(1);%mmHg
VT0init=ParamICvals(2);%ml
ChSlinit=ParamICvals(3);%mmHg/ml
XJvinit=ParamICvals(4);% nd
XJsinit=ParamICvals(5);% nd
Ppinit=ParamICvals(6);% mmHg
phimaxinit=ParamICvals(7);%ml/min
PosSlinit=ParamICvals(8);%mmHg

baselineValues=[P0init VT0init ChSlinit XJvinit XJsinit Ppinit ...
    phimaxinit PosSlinit];% parameter default initial guesses

%parameter and initial conditions for simulated results
% Possenti parameters
p50=2;%mmHg
piinterval=4-(-2);% range of interstitial pressures for sigmoidal variation
phimin=2;% scaled min lymph flux

% baseline parameters; perhaps perturbed from the values set above
P0=baselineValues(1);%interstitial pressure breakpoint from Chapple 1993 Fig A2
VT0=baselineValues(2);%Chapple breakpoint volume
ChSl=baselineValues(3);%Chapple slope
XJv=baselineValues(4);
XJs=baselineValues(5);
Pp=baselineValues(6);
phimax=baselineValues(7);
PosSl=baselineValues(8);

pT0init=chapplepressure(Vi0);%Chapple 1993 Fig A2

IC=[Vp0 Cp0 Vi0 Ci0];
H0 = (1-Hct0)/(Hct0*Vp0); %constant for Hct (t) calculation

ap1 = 0.1752; %coef for Cp to get pi_p in mmHg(ml/mg)
ap2 = 0.0028; %coef for Cp^2 
ai1 = 0.2336; %coef for Ci to get pi_i in mmHg(ml/mg)
ai2 = 0.0034; %coef for ci^2 

pip=(ap1+ap2*Cp0)*Cp0;% oncotic pressures
pii=(ai1+ai2*Ci0)*Ci0;

% Definition of constantset--contains some archaic, unused values
constantset(1)=ap1;constantset(2)=ap2;constantset(3)=ai1;constantset(4)=ai2;
constantset(7)=Hct0;constantset(8)=P0;constantset(9)=ChSl;constantset(10)=VT0;
constantset(11)=phimax;constantset(12)=phimin;constantset(13)=PosSl;
constantset(14)=p50;

pressures=[Ppinit pT0init pip pii];

% These are the regression coefficients obtained for set values of the
% permeability regression coefficients
JvNcoefs=[3.11850348299278;58.2668067303779;-58.2684565822765;...
    -47.3109808554083;40.7165353119554];
JsNcoefs=[-2917.94139376802;131.869569261009;-131.920965019932;...
    96.1472683480481;49.6265280379148];

% % set bounds for parameter search range
% %ranges of acceptable starting points and [baseline values]
LoSet=LoBnds;
HiSet=HiBnds;
lb=LoSet(paramidsToOpt);
ub=HiSet(paramidsToOpt);
%% Construct sets of random starting points for optimization

pstart=zeros(length(nparams2opt),Ntrials);
for j=1:nparams2opt% construct LHS samples of random starting points
    pstart(j,:)=LHS_Call(lb(j),baselineValues(j),ub(j),0,Ntrials,'unif');
end

disp(approachinfo);
disp(ModelID);
disp(datafilename);
disp(solverinfo);
disp(allParamNames);

%% Now do optimization
solveropts=odeset('RelTol',1e-6,'MaxStep',deltat);% set options for ode solver
% options for search process
searchOptions=optimset('Display','notify','TolFun',1e-9,...
    'TolX',1e-9,'MaxIter',5000,'MaxFunEvals',10000);

%initialize all histories though only some will be used
P0history=zeros(Ntrials,1);% preallocate results storage vectors
XJvhistory=zeros(Ntrials,1);
XJshistory=zeros(Ntrials,1);
Pphistory=zeros(Ntrials,1);
Cphistory=zeros(Ntrials,1);
Cihistory=zeros(Ntrials,1);
varhistory=zeros(Ntrials,1);

tic
starttime=tic; %start timer
parfor i=1:Ntrials %parallel optimization from each initial starting set
    initialparamguess=[];
    for j=1:length(paramidsToOpt) %set as many initial conditions as needed
        initialparamguess=[initialparamguess pstart(j,i)];
    end
    %initialparamguess
    x=initialparamguess;%initial values for optimized parameters

    %construct melded list of parameter values 
    allParamValues=baselineValues;% this holds initial values for all parameters
    allParamValues(paramidsToOpt)=pstart(:,i);%these are initial values for melded set

    f=@(x) ssq_revstarlingoptRRI(x,tfit,hctfit,...
        Hct0,paramidsToOpt,allParamValues,tspanfit,IC,...
        solveropts,interpApplFlux,constantset);

    [x,fval,exitflag,output] = fminsearchbnd(f,initialparamguess,lb,ub,searchOptions);
    x;% vector of optimized parameter values
    output;

    for j=1:length(x) % update the histories of the optimized parameters
        switch paramidsToOpt(j)
            case 1
                P0history(i)=x(j);% save optimal results
            case 2
                VT0history(i)=x(j);
            case 3
                ChSlhistory(i)=x(j);
            case 4
                XJvhistory(i)=x(j);
            case 5
                XJshistory(i)=x(j);
            case 6
                Pphistory(i)=x(j);
            case 7
                phimaxhistory(i)=x(j);
            otherwise
                PosSlhistory(i)=x(j);
        end
    end

    varhistory(i)=fval;

end% of i=1:Ntrials loop

[M,I]=min(varhistory);% find parameter set I corresponding to best fit

% do a final search using the min variance parameter value set

nearlyOptimizedParamValues=[];
for j=1:length(paramidsToOpt)
    switch paramidsToOpt(j)
        case 1
            nearlyOptimizedParamValues(j)=P0history(I);% optimal P0
        case 2
            nearlyOptimizedParamValues(j)=VTOhistory(I);
        case 3
            nearlyOptimizedParamValues(j)=ChSlhistory(I);
        case 4
            nearlyOptimizedParamValues(j)=XJvhistory(I);
        case 5
            nearlyOptimizedParamValues(j)=XJshistory(I);
        case 6
            nearlyOptimizedParamValues(j)=Pphistory(I);
        case 7
            nearlyOptimizedParamValues(j)=phimaxhistory(I);
        otherwise
            nearlyOptimizedParamValues(j)=PosSlhistory(I);
    end
end


%merge nearly optimized parameter values with rest of parameter values
nearlyFinalParamValues=baselineValues;% this holds initial values for all parameters
nearlyFinalParamValues(paramidsToOpt)=nearlyOptimizedParamValues;%merge optimized values

%Redo optimal search from minimum variance trial in set
f=@(x) ssq_revstarlingoptRRI(x,tfit,hctfit,...
    Hct0,paramidsToOpt,nearlyFinalParamValues,tspanfit,IC,...
    solveropts,interpApplFlux,constantset);

[x,fval,exitflag,output] = fminsearchbnd(f,nearlyOptimizedParamValues,lb,ub,searchOptions);

%% Now compute fit

optimizedParamValues=x;% best values after final search

%merge optimized parameter values with rest of parameter values
finalParamValues=baselineValues;% this holds initial values for all parameters
finalParamValues(paramidsToOpt)=optimizedParamValues;%merge optimized values

P0opt=finalParamValues(1);VT0opt=finalParamValues(2);ChSlopt=finalParamValues(3);
XJvopt=finalParamValues(4);XJsopt=finalParamValues(5);Ppopt=finalParamValues(6);
phimaxopt=finalParamValues(7);PosSlopt=finalParamValues(8);

tspanuse=tspanall;% use optimal parameters over entire range of data
%IC=[Vp0 Cp0 Vi0 Ci0];

sol = ode15s(@(t,y) revstarlingmodelRRI(t,y,interpApplFlux,constantset,...
    finalParamValues),tspanuse,IC,solveropts);

Y = deval(tall,sol);%gets vectors of state variables at tall times
Pvol = Y(1,:);%plasma volume
Vpbest=Pvol;
Vp0=Vpbest(1);
Cpbest=Y(2,:);%plasma protein concentration
Vibest=Y(3,:);%interstitial volume
Cibest=Y(4,:);%interstitial protein concentration
H0 = (1-Hct0)/(Hct0*Vp0); %constant for Hct (t) calculation
Hct = 1./(H0*Pvol+1); % fitted time course of Hct using minimum variance values

finalvar=var(Hct-hct);
formatspec3='Ending minimum variance = %8.6e \n';
fprintf(formatspec3,finalvar);

%% compute fluxes for optimized model

P0final=finalParamValues(1);VT0final=finalParamValues(2);
ChSlfinal=finalParamValues(3);
XJvfinal=finalParamValues(4);XJsfinal=finalParamValues(5);
Ppfinal=finalParamValues(6);phimaxfinal=finalParamValues(7);
PosSlfinal=finalParamValues(8);

%reconstruct min variance vectors of oncotic pressures
pip_final = ap1*Cpbest + ap2*Cpbest.^2; %Eqn 3
pii_final = ai1*Cibest + ai2*Cibest.^2; %Eqn 3

%preallocate vector space
Jv=zeros(1,length(pip_final));%vector for piP(t)
Js=zeros(1,length(pii_final));%vector for piT(t)
Jl=zeros(1,length(pii_final));%vector for Jl

for i=1:length(t) %find Jv and Js from [pP,pT,piP, piT] (-:out of plasma)
    pTifinal(i)=chapplepressure(Vibest(i));%interstitial hydrostatic pressure
    pressuresi=[Ppfinal pTifinal(i) pip_final(i) pii_final(i)];%boundary conditions for ith point
    % the (-) for Jv and Js converts Facchini to to delosReyes convention,
    % the latter being flow INTO the plasma is positive
    Jv(i)=-XJvfinal*(JvNcoefs(1)+pressuresi*JvNcoefs(2:5));
    Js(i)=-XJsfinal*(JsNcoefs(1)+pressuresi*JsNcoefs(2:5));
    % find Jl from Possenti using Chapple
    denom1=1+exp((pTifinal(i)-p50)/(PosSlfinal));
    Jl(i)=phimaxfinal-((phimaxfinal-phimin)/denom1);%Possenti eqn 4, pg 103
end

%% Output results
% outputs:
% t= time vector, same as input,
% FittedHct= fitted hct(t) using optimized parameter values,
% FittedSVs= fitted state variables [Vp Cp Vi Ci],
% FittedParamValues= optimized parameter values for the fitted parameters,
% FittedPressures= hydrostatic and oncotic pressure for optimal model,
% FittedFluxes= Jv, Js and Jl for optimal model,
% TimeDate= time and date stamp when optimization was done, and
% FittedVar= variance of the fit

% % Foutputs: 

Fresults.t=t;
Fresults.fittedHct=Hct;
Fresults.fittedSVs=[Vpbest;Cpbest;Vibest;Cibest];
Fresults.fittedParamValues=finalParamValues;
Fresults.fittedPressures={Pp,pTifinal,pip_final,pii_final};%cell due to Pp
Fresults.fittedFluxes=[Jv;Js;Jl];
Fresults.TimeDate=TimeDate;
Fresults.MinVar=finalvar;

optResultsOutputRRI %plots various results

if fileSaveFlag==1 %save results to file
    optResultsSaveOutputRRI
end

toc
%--------------------------------------------
end