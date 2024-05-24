function [ssqdevs] = ssq_revstarlingoptRRI(x,tfit,hctfit,...
        Hct0,paramids2opt,allParamValues,tspanuse,IC,...
        solveropts,interpApplFlux,constantset) 

%  computes ssq deviations between model and data to use as the objective
%  function for use with fminsearch parameter estimation.
%  3/16/24 adapted to use net applied flux into plasma

% x is the vector of parameter values; e.g., x = [XJv Xjs Pp] 
params2opt=x;% THIS IS A CRITICAL STATEMENT to force optimization to proceed

ap1=constantset(1);ap2=constantset(2);ai1=constantset(3);ai2=constantset(4);

paramstoUse=params2opt;%this will hold parameter values used in search
j=1;%index for paramsinSet
for i=1:length(allParamValues) % construct the meld of optimized and fixed parameters
    if j>length(paramids2opt)% there are no more parameters to update
        paramstoUse(i)=allParamValues(i);%so fill out the parameter list
    elseif i==paramids2opt(j) %j<=#parameters to optimize
        paramstoUse(i)=x(j);
        j=j+1;
    else %j<=# parameters to update but there are more to add to the list
        paramstoUse(i)=allParamValues(i);
    end
end
paramset=paramstoUse;

ssqdevs=0; %initialize ssd

% find solution over times tspanuse: fit times or all times

sol = ode15s(@(t,y) revstarlingmodelRRI(t,y,interpApplFlux,constantset,paramset),tspanuse,IC,solveropts);

% Protect against solution terminating due to failure to meet tolerance
% by truncating tfit to the highest time point available in sol.x
isolx=find(tfit<=sol.x(end),1,'last'); %highest index in tfit that can be used
if isolx~=length(tfit)
    isolx 
end
tfitsolx=tfit(1:isolx);% truncates time for deval
Y = deval(tfitsolx,sol);%gets vectors of state variables at tfit times
 
Vp = Y(1,:);%plasma volume for valid time span

% Vp0=Vp(1);
H0 = (1-Hct0)/(Hct0*Vp(1)); %update H0 constant for Hct (t) calculation
%   from new plasma volume; this allows Vp to be optimized
Hct = 1./(H0*Vp+1); % time course of Hct over valid time span
hctfitsolx=hctfit(1:isolx);% hct over valid tome span
if size(Hct)~=size(hctfitsolx) %ensure Hct is a row vector to match hctfit
    Hct=Hct';
end
devs=(Hct-hctfitsolx);%vector of differences between model and observations

ssqdevs=sum(devs.^2);