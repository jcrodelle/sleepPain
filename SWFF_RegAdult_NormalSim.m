function y = SWFF_RegAdult_NormalSim(t, x)

% SWFF model with parameters based on Piltz, Diniz Behn, Booth JTB 2020
% Regular sleepers, median parameters

% time scale in minutes
fW = x(1);
fS = x(2);
fSCN = x(3);
h = x(4);
% for the oscillator model for circadian clock
c = x(5);
xc = x(6);
n = x(7);
  

% SW Parameters, regular median values
gsw=0.2508;
gscnw= 0.01;  
gws=0.25;
gscns= 0.07; 

% parameter values from JTB Piltz paper -- human "regular sleeper" data 
tauW=23;
tauS=10;
maxW=6;
maxS=6;
betaW=-0.4;
alphaW=0.4;
alphaS=0.2;
k1 = -0.0118;  
k2 = -0.005;       
thetaW=4;
tauhw=946.8; 
tauhs=202.2;
hmax=323.88;
hmin=0;

maxSCN=7;
betaSCN=-0.1;
alphaSCN=0.7;
tauSCN=0.5;


% Clock Parameters
 taux=24.2;
 alph_0=0.05;
 beta=0.0075;
 G=33.75;
 p=0.5;
 k=0.55;
 I_0=9500;
 mu=0.23;
 
% homeostatically varying S threshold
betaS=k2*h+k1;
 

% % -------------------------------------------------
% 
% %circadian drive from Forger, Jewett, Kronauer 1999 model
% 
% % Light function I 
% 
% % I as constant

% % I fluctuating repeating step function: L:D = 16:8
% t = 0 at beginning of light period (8 am)
per=24*60;% minutes -- 24 hours
t2=12*60; % Light duration duration 
Iamp=600; % Lux value -- 


I=Iamp*heaviside(-1*((t-floor(t/per)*per)-t2));


% inputs to circadian clock model
% %Drive due to Process L  
alph=(alph_0)*((I^p)/((I_0)^p));
B_hat=G*(1-n)*alph;

% Sensitivity modulation
B=B_hat*(1-0.4*c)*(1-0.4*xc);


%inputs for steady state response functions
W_INPUT = gscnw*fSCN - gsw*fS;
SCN_INPUT = c;
S_INPUT = -gws*fW - gscns*fSCN;

%t>t_sd1 & t < t_sd1+length_sd

%population steady state response functions
% add 10 here
W_inf = maxW*0.5*(1+tanh((W_INPUT-betaW)/alphaW));
SCN_inf = maxSCN*0.5*(1+tanh((SCN_INPUT-betaSCN)/alphaSCN));
S_inf = maxS*0.5*(1+tanh((S_INPUT-betaS)/alphaS));

%ODES

y = zeros(7,1);

y(1) = (W_inf-fW)/tauW; 
y(2)=(S_inf-fS)/tauS;
y(3)=(SCN_inf-fSCN)/tauSCN;
y(4)=heaviside(fW-thetaW)*(hmax-h)/tauhw+heaviside(thetaW-fW)*(hmin-h)/tauhs;

% Clock equations

y(5)=(pi/720)*(xc+B);
y(6)=(pi/720)*(mu*(xc-((4*xc^3)/3))-c*(((24/(0.99669*taux))^2)+k*B));
y(7)=((alph*(1-n))-(beta*n));


end
