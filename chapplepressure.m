function [Pi]=chapplepressure(Vi)
%returns interstitial pressure as a function of interstitial volume
%according to data from Chapple Fig.A.2 from Chapple et al: A model
% of human microvascular exchange: parameter estimation based on normals 
% and nephrotics;Computer Methods and Programs in Biomedicine 41 (1993) 33-54
% Vi in ml, Pi in mmHg
% 2/9/23

if Vi<9500 %stiff compliance below Vi breakpoint, normal volemic
    P0=-7.2; V0=5000;K=1.82e-3;% K is slope in mmHg/ml
elseif Vi>=9500 && Vi<10300 %transition Vi btwn stiff and compliant regions
    P0=1; V0=9500;K=1.25e-3;
else
    P0=2;V0=10300;K=5e-5;% hypervolemic region
end

Pi=P0+(Vi-V0)*K;
end
