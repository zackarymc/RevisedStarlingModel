function [JvN,JsN]=facchinifluxcalc(pressures)
% Computes JvN and JsN at a single point in time for given pressures 
% values using Facchini regression coefficients and Possenti parameters.
% Use this in a loop to generate a time vector of fluxes.

% This differs from facchinifluxcalxolder.m by having the regression
% coefficients embedded in the function, so the only calling parameters
% required are the pressures = [pP,pT,piP,piT]

% JvN and JsN are the normalized fluxes of water and proteins
% respectively at the point in time t(i)

%find JvN and JsN from [pP,pT,piP, piT] (-:out of plasma)
% pressures=[pG0best pTi(i) pi_pbest(i) pi_ibest(i)]=boundary conditions
% for ith point.
% This routine returns normalized fluxes in the Facchini convention in  
% which flow INTO the plasma is positive

% NOTE: JvN > 0 and JsN > 0 are computed as the fluxes  out of the plasma 
% and into the interstitium.
% The de los Reyes convention mandates a (-) conversion for Jv and Js as 
% JvN and JsN are used later in the model.

% Facchini regresion constants
% These are the regression coefficients obtained for set values of the
% permeability regression coefficients
JvNcoefs=[3.11850348299278;58.2668067303779;-58.2684565822765;...
    -47.3109808554083;40.7165353119554];
JsNcoefs=[-2917.94139376802;131.869569261009;-131.920965019932;...
    96.1472683480481;49.6265280379148];

JvN=(JvNcoefs(1)+pressures*JvNcoefs(2:5));% water flux out of the plasma
JsN=(JsNcoefs(1)+pressures*JsNcoefs(2:5));%protein flux out of the plasma

end

