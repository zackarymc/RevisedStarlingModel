function ydot=revstarlingmodelRRI(t,y,interpApplFlux,constantset,...
paramset)
% This solves the system described in the paper by A.A de los Reyes V and
% others with Peter Kotanko with added glycocalyx factor using the results
% of Facchini et al 2014 implemented using regression results. Possenti
% computation of lymphatic flow is used with Pi computed ala Chapple 1993.

% Calling arguments are t: time for solution, y:state variables, 
% interpApplFlux:net flux into plasma at solution times,
% newconstantset: coefs to compute oncotic pressures from
% concentrations, paramset: all parameter values to use in solution.

% Lymphatic flow (Jl) calculated using Possenti et al., Numerical 
% simulations of the microvascular fluid balance with a non-linear model 
% of the lymphatic system. Microvascular Research 122:101-110, 2019.
% This approach uses a sigmoidal function shown in Fig. 2 with Eq. 4 with
% parameters listed on page 103.

% Interstitial pressure calculated using Chapple et al., Computer Methods
% and Programs in Biomedicine 41: 33-54, 1993 using 3 segment linear fits 
% to data from Fig. A2, pg. 52 for Vi<9.5, 9.5<Vi<10.3, and Vi>10.3 L.

% 1/30/2024: cleaned up code and comments

% Generate JApplFlux from at solver time t
JApplFlux=interpApplFlux(t);

%The state and auxiliary equations are coded below for the four state
%variables: y(1) = Vp (volume of plasma compartment)
%           y(2) = Cp (concentration of protein in plasma comapartment)
%           y(3) = Vi (volume of interstitiall compartment)
%           y(4) = Ci (concentration of protein in interstitial compartment)

Vp=y(1);% state variables
Cp=y(2);
Vi=y(3);
Ci=y(4);

P0=paramset(1);VT0=paramset(2);ChSl=paramset(3);XJv=paramset(4);
XJs=paramset(5);pP=paramset(6);phimax=paramset(7);
PosSl=paramset(8);

% Definition of constantset
% constantset(1)=ap1; constantset(2)=ap2; constantset(3)=ai1;
% constantset(4)=ai2;
ap1=constantset(1);ap2=constantset(2);ai1=constantset(3);
ai2=constantset(4);

%reconstruct plasma and interstitial oncotic pressures from concentrations
piP = ap1*Cp + ap2*Cp.^2; %Eqn 3
piT = ai1*Ci + ai2*Ci.^2; %Eqn 3

% Interstitial pressure estimated from Chapple 1993, Fig A2 and text
pT=chapplepressure(Vi);%interstitial hydrostatic pressure, smoothed Chapple

% Possenti parameters
phimin=2;
p50=2;

% boundary conditions for Facchini model regression
pressures=[pP pT piP piT];

[JvN,JsN]=facchinifluxcalc(pressures);% These are the normalized water and 
% protein fluxes according to Facchini's convention in which flows into the  
% interstitial space out of the plasma are positive and due solely to the 
% pressures (hydrostatic and oncotic) working on the glycocalyx.
% Any lymphatic flux into plasma is additional. 
% Ergo: Facchini flux convention is fluxes>0 into Vi and out of Vp.

% Convert Facchini fluxes to de los Reyes fluxes where fluxes>0->Vp.
Jv=-XJv*JvN;% actual water flux in dlR convention where plasma inflow>0
Js=-XJs*JsN;% JsN>0 is into Vp (Facchini) so Js<0 is out of Vp (dlR)
% Calculate Jl ala Possenti eqn 4, pg 103
denom1=1+exp((pT-p50)/(PosSl));
Jl=phimax-((phimax-phimin)/denom1);%Possenti 

Jsnet=Js+Jl*Ci;%total protein flux Vi->>Vp (Jl>0 and Js>0 go into plasma)

% The following state equations are from Reyes et al., Eqn 9 and follow the
% dlR convention of fluxes > 0 are INTO the plasma
dVpdt = Jv + Jl + JApplFlux;% total water flux into plasma
dCpdt = (Jsnet - Cp*(dVpdt))/Vp;% plasma protein
dVidt = -Jv - Jl;% interstitial volume flux balance
dCidt = (-Jsnet + Ci*(-dVidt))/Vi;% interstitial protein flux balance
ydot(1) = dVpdt; 
ydot(2) = dCpdt;
ydot(3) = dVidt;
ydot(4) = dCidt;
ydot=ydot'; %get column vector for ode solver

end
%----end of function revstarlingmodelRRI-------------------
