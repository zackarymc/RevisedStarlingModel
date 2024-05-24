% script optResultsOutputRRI.m to plot results of selected optimization

% Instigated from function revStarOptSelectedRRI.m
numfig=1;%initialize figure numbering
%% Plot results
% Results vs data: "FIT"

% plot results from best minimum variance results after final search versus reference data

figure(numfig);clf % Plot reference data used for the optimization and checking
numfig=numfig+1;% figure number for next figure

yyaxis left
plot(tfit,hctfit,'-k','LineWidth',.2) %plots data results
hold on
plot(tcheck,hctcheck,'-k','LineWidth',.2) %plots
hold on
xlabel('time (min)');
ylabel('Hct');hold on;

plot(tall,Hct,':c','LineWidth',4);hold on %plots final fitted results

yyaxis right
plot(t,appliedFlux,'k--','LineWidth',3)
ylabel('Applied Flux (ml/min)')

legend('Hct','fit','Applied Flux (ml/min)')

formatSpec="Fit after final fminsearch for %s \n on %s for %s with %i trials, using %i datapoints";
A2=approachinfo;A3=ModelID;A4=char(datetime);A5=patientID;A6=ndatatotal;A1=Ntrials;
str=sprintf(formatSpec,A3,A4,A5,A1,A6);
title(str);

HctEst0=Hct(1);% calculate initial and last Hct and volumes and deltas
HctEstlast=Hct(end);
PvolEst0=Pvol(1);
PvolEstlast=Pvol(end);
Ivol=Y(3,:);%interstitial volume
IvolEst0=Ivol(1);
IvolEstlast=Ivol(end);
TotVolEst0=PvolEst0+IvolEst0;
TotVolEstlast=PvolEstlast+IvolEstlast;
dPvol=PvolEstlast-PvolEst0;
dIvol=IvolEstlast-IvolEst0;
dTotVol=TotVolEstlast-TotVolEst0;

dHct=hctall-Hct;%deviations btwn data and final search fit
bestfittedvariance=var(dHct);% best formulation for variance;
finalfitvar=num2str(round(bestfittedvariance,9));%from final search

optvarnames=char(allParamNames(paramidsToOpt));
optvarstring=['var'];

for i=1:length(paramidsToOpt)
    optvarstring=[optvarstring,'/',strtrim(optvarnames(i,:))];
end

s1=['Optimized ',optvarstring,'=',finalfitvar,'/'];

for i=1:length(paramidsToOpt)
    s1=[s1,num2str(optimizedParamValues(i)),'/'];
end

dim = [.15 .6 .3 .3];
annotation('textbox',dim,'String',s1,'FitBoxToText','on');
figure(numfig-1);%this is somehow needed to be able to save the correct figure???

%% plot 2x2 summary plot

figure(numfig);clf
numfig=numfig+1;
subplot(2,2,1)
yyaxis left
plot(tall,Vpbest,'Linewidth',2)
xlabel('Time in min')
ylabel('Plasma Vol in ml')
title('Actual Responses')
yyaxis right
plot(tall,appliedFlux,'--','Linewidth',2)
ylabel('Applied Flux in ml/min')
subplot(2,2,2)
yyaxis left
plot(tall,Cpbest,'Linewidth',2)
ylabel('Plasma Protein conc in mg/ml')
xlabel('Time in min')
yyaxis right
plot(tall,Cibest,'--','Linewidth',2)
ylabel('Interstitial Protein conc in mg/ml')
%title('Interstitial protein')
subplot(2,2,3)
plot(tall,Vibest,'Linewidth',2)
ylabel('Interstitial volume in ml')
xlabel('Time in min')
%title('interstitial volume')
subplot(2,2,4)
plot(tall,hctall,'k--');hold on
plot(tall,Hct,'r','LineWidth',4);hold on
xlabel('Time in min')
ylabel('Hct')
%title('Time course of Hct')

s1=['variance = ',num2str(bestfittedvariance)];
dim = [.6 .15 .3 .3];
annotation('textbox',dim,'String',s1,'FitBoxToText','on');

s2=[];
for i=1:length(paramidsToOpt)-1
    s2=[s2,num2str(optimizedParamValues(i),4),'/'];
end
s2=[s2,num2str(optimizedParamValues(i+1),4)];

dim2=[0.139285714285714,0.103571428571429,0.34375,0.053571428571429];
annotation('textbox',dim2,'String',s2,'FitBoxToText','on');
%% Now plot histograms of parameter  starting points and final fits
if Ntrials>20 % don't plot if too few trials
    figure(numfig);clf % Plot starting point histograms
    numfig=numfig+1;% figure number for next figure

    for i=1:nparams2opt % do histograms for optimized parameters only
        subplot(nparams2opt,1,i)
        histogram(pstart(i,:))
        if i==1
            title(['Distribution of starting values for ',num2str(Ntrials),' trials'])
        end
        ylabel(optvarnames(i,:))
    end

    figure(numfig);clf % Plot parameter histograms
    numfig=numfig+1;% figure number for next figure
    optparamhistory=zeros(nparams2opt,Ntrials);

    for j=1:nparams2opt% select the opt values to plot
        switch paramidsToOpt(j)
            case 1
                optparamhistory(j,:)=P0history;% optimal P0
            case 2
                optparamhistory(j,:)=VT0history;
            case 3
                optparamhistory(j,:)=ChSlhistory;
            case 4
                optparamhistory(j,:)=XJvhistory;
            case 5
                optparamhistory(j,:)=XJshistory;
            case 6
                optparamhistory(j,:)=Pphistory;
            case 7
                optparamhistory(j,:)=phimaxhistory;
            otherwise
                optparamhistory(j,:)=PosSlhistory;
        end
    end

    for i=1:nparams2opt % do histograms for optimized parameters only
        subplot(nparams2opt,1,i)
        lowbin=min(optparamhistory(i,:));hibin=max(optparamhistory(i,:));
        %histogram(optparamhistory(i),'BinEdges',[lowbin,hibin])
        histogram(optparamhistory(i,:),'NumBins',20)
        if i==1
            title(['Distribution of optimized values for ',num2str(Ntrials),' trials'])
        end
        ylabel(optvarnames(i,:))
    end
end
%--------------------------------------------
%%
%plot time courses of Vp/Cp/Vi/Ci/Juf for best parameters
figure(numfig);clf % Plot time courses of state variables and forcing function
numfig=numfig+1;% figure number for next figure

subplot(5,1,1) %plot Vp
plot(tall,Vpbest);hold on
title('State variable time courses for minimum variance estimates')
xlabel('time (min)')
ylabel('Plasma volume (ml)')

subplot(5,1,2)
plot(tall,Cpbest);hold on
xlabel('time (min)')
ylabel('Plasma protein conc (mg/ml)')

subplot(5,1,3)
plot(tall,Vibest);hold on
xlabel('time (min)')
ylabel('Interstitial volume (ml)')

subplot(5,1,4)
plot(tall,Cibest);hold on
xlabel('time (min)')
ylabel('Interstitial protein conc (mg/ml)')

subplot(5,1,5)
plot(tall,appliedFlux)
xlabel('time (min)')
ylabel('Applied Flux (ml/min)')

%% Plot fluxes for optimized mpdel
% time courses of Jv/Js/Jl/Juf for optimized parameters
figure(numfig);clf % Plot 3D space of optimal values
numfig=numfig+1;% figure number for next figure

subplot(5,1,1)
plot(tall,Jv)
xlabel('time (min)')
ylabel('Jv (ml/min)')
title('Fluxes')

subplot(5,1,2)
plot(tall,Js)
xlabel('time (min)')
ylabel('Js (mg/min)')

subplot(5,1,3)
plot(tall,appliedFlux)
xlabel('time (min)')
ylabel('Applied Flux (ml/min)')

subplot(5,1,4)
plot(tall,Jl)
xlabel('time (min)')
ylabel('Jl (ml/min)')

subplot(5,1,5)
plot(tall,pTifinal)
xlabel('time (min)')
ylabel('pT (mmHg)')

%% do flux check if desired
if fluxcheckflag==1 % do separate check of volume changes to verify fluxes
    %find sum of fluxes (i.e. net volume change) in ml into plasma or into
    %interstitium

    %integrate fluxes to get net volume changes
    summedJv=simps(t,Jv)% Jv>0 into plasma
    summedJl=simps(t,Jl)% Jl>0 into plasma
    summedAppliedFlux=simps(t,appliedFlux)% Juf<0 out of plasma
    summedJs=simps(t,Js)% Js>0 into plasma
    summedJlxCi=simps(t,Jl.*Cibest)% JlxCi>0 into plasma always the case
    summedJinf=simps(t,Jinf)
    summedJuf=simps(t,Juf)

    summedJvin=0;summedJvout=0;summedJsin=0;summedJsout=0;
    summedAppliedFluxin=0;summedAppliedFluxout=0;
    summedJlin=0; summedJlout=0;

    for i=1:length(t)
        if Jv(i)>0
            summedJvin=summedJvin+Jv(i);%fluid into plasma
        else
            summedJvout=summedJvout+Jv(i);%fluid out of plasma into interstitium
        end
        if Js(i)>0
            summedJsin=summedJsin+Js(i);%protein into plasma
        else
            summedJsout=summedJsout+Js(i);%protein out of plasma into interstitium
        end
        if appliedFlux(i)>0;%net flux into plasma
            summedAppliedFluxin=summedAppliedFluxin+appliedFlux(i);
        else
            summedAppliedFluxout=summedAppliedFluxout+appliedFlux(i);
        end
        if Jl(i)>0%lymph flow into plasma
            summedJlin=summedJlin+Jl(i);
        else
            summedJlout=summedJlout+Jl(i);
        end
    end

    % compute separately the volumes of fluid and protein
    % plasma fluid
    summedJvin=summedJvin*deltat%fluid into plasma ml
    summedJvout=summedJvout*deltat%fluid out of plasma into interstitium ml

    % plasma protein
    summedJsin=summedJsin*deltat%protein into plasma mg
    summedJsout=summedJsout*deltat%protein out of plasma into interstitium mg

    % net applied fluxes
    summedAppliedFluxin=summedAppliedFluxin*deltat
    summedAppliedFluxout=summedAppliedFluxout*deltat

    %net lymphatic fluxes
    summedJlin=summedJlin*deltat;
    summedJlout=summedJlout*deltat;

    Vpend=Vpbest(end);Vpstart=Vpbest(1);Cpend=Cpbest(end);Cpstart=Cpbest(1);
    Viend=Vibest(end);Vistart=Vibest(1);Ciend=Cibest(end);Cistart=Cibest(1);

    % compute the changes in plasma volumes
    % fluid volume changes
    netdeltaVpfluxes=summedJv+summedJl+summedAppliedFlux%change in plasma fluid volume (ml)
    dVpvol=Vpend-Vpstart%actual Vp change from start
    Vpdifference_Fluxes_minus_dVpvol=netdeltaVpfluxes-dVpvol

    % compute the changes in interstitial volumes
    % fluid volume changes
    netdeltaVifluxes=-summedJv-summedJl%change in interstitial fluid volume (ml)
    dVivol=Viend-Vistart%actual Vi change from start
    Vidifference_Fluxes_minus_dVivol=netdeltaVifluxes-dVivol

    % plasma protein changes
    % summed protein flux into PLASMA
    % NOTE: dlR accounts for lymphatic flux in calculation for Js
    Net_lymphatic_protein_flux=summedJlxCi%protein flux via lymph flow out of Vi into Vp
    net_plasma_protein_change_fluxes=Net_lymphatic_protein_flux + summedJs%total protein out of Vi
    dPlasmaProtein_volXconc=Vpend*Cpend-Vpstart*Cpstart
    PlasmaProteinDifference=net_plasma_protein_change_fluxes-dPlasmaProtein_volXconc

    % interstitial protein changes
    % summed interstitial fluxes into interstitium
    % net_interstitial_protein_fluxes=sum(fluxes into Vi) so -summedJs
    %   since Js<0 is flux into Vi, and Jl*Ci is out of Vi
    net_interstitial_protein_fluxes=-summedJs-Net_lymphatic_protein_flux
    % net_interstitial_protein_fluxes should = ∆Vol*Conc, end - start
    dInterstitialProtein_volXconc=Viend*Ciend-Vistart*Cistart
    % But the above calculations have opposite signs!!! SO the next line is
    %   2x instead of ~ 0.
    InterstitialProteinDifference=net_interstitial_protein_fluxes-dInterstitialProtein_volXconc

    % Plasma protein amounts
    PlasmaProteinBefore=Cpstart*Vpstart
    PlasmaProteinAfter=Cpend*Vpend
    deltaPlasmaProtein=PlasmaProteinAfter-PlasmaProteinBefore

    %interstitial protein amounts
    InterstitialProteinBefore=Cistart*Vistart
    InterstitialProteinAfter=Ciend*Viend
    deltaInterstitialProtein=InterstitialProteinAfter-InterstitialProteinBefore

    % Total protein amounts
    TotalMass_before=Cpstart*Vpstart+Cistart*Vistart
    TotalMass_after=Cpend*Vpend+Ciend*Viend
    TotaldeltaProteinMass=TotalMass_after-TotalMass_before

    % Volume conservation
    dVplasma=Vpend-Vpstart
    summedAppliedFlux
    dVinterstitium=Viend-Vistart
end

optvarstring %print optimized parameter names
finalfitvar
"final parameter values:"

