% RRIreconstruct.m
% script to construct evenly sampled hct(t) and t data plus juf(t) and
% Jinf(t) for RRI data from graphical digitizations

% IMPORTANT!
% Use the Import procedure to input t and hct data from WebPlotDigitizer 
% Excel file previously saved. This gives an nx2 table with t in column1
% and hct in column 2.

% Before running, be sure to set the individual Jinf parameters in the
% below location

% ONLY THEN RUN THE SCRIPT
% executable code follows:

S=whos;
fName=S.name;%recovers filename for imported table
tvarname=string(strcat(string(fName),'.t'));% name of tabledata for t
torig=eval(tvarname);% gets t data
hctvarname=string(strcat(string(fName),'.hct'));% name of tabledata for hct
hctorig=eval(hctvarname);%gets hct data

% Set these parameters for each individual data set
%  for zdv250960
% jufAmplitude=14.7;% ml/min
% infusionTime=[0 5.51];% min
% infusionAmplitude=[250 250]; % ml 
% infusionDuration=[5 5];% min
% infusionRate=infusionAmplitude./infusionDuration;

%  for zdv2501058
jufAmplitude=8.5;% ml/min
infusionTime=[0 6.01];% min
infusionAmplitude=[250 250]; % ml 
infusionDuration=[5 5];% min
infusionRate=infusionAmplitude./infusionDuration;

dt=1;%sampling interval (min)
sf=fix(1/dt);

t=0:dt:torig(end);
hct=interp1(torig,hctorig,t,'spline');
N=length(t);

Juf=jufAmplitude*ones(1,N);%constant juf
Jinf=zeros(1,N);
for i=1:N% create Jinf
    if (t(i)>=infusionTime(1)) && (t(i)<(infusionTime(1)+infusionDuration(1)))
        Jinf(i)=infusionRate(1);
    elseif (t(i)>=infusionTime(2)) && (t(i)<(infusionTime(2)+infusionDuration(2)))
        Jinf(i)=infusionRate(2);
    else
        Jinf(i)=0;
    end
end

Japplied=Jinf-Juf;% net flux into/out of plasma

figure(1);clf
yyaxis left
plot(torig,hctorig,'o')
hold on
plot(t,hct,'-x')
xlabel('time (min)')
ylabel('hct')

yyaxis right
plot(t,Japplied)
ylabel('net applied flux')

newfName=strcat(fName,num2str(sf));

save(newfName,'t','hct','Juf','Jinf')

