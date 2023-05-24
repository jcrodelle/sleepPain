% This script will run a normal day and compute sleep and wake times
% from "Modeling homeostatic and circadian modulation of human pain
% sensitivity" Crodelle, Vanty, Booth 2023
clear all
close all

d = 7; %number of days
min_d = d*24*60; % in minutes
dt = 0.25; % time step in minutes
tSpan=0:dt:min_d; % span of time for simulation
% initial conditions
ics =[4.0865    0.0657    2.9271   26.4815   -0.2130    1.0356    0.5755];
% for determining the timing of sleep and wake
opts = odeset('Events', @events);
% run the ODE solver
[t,x,te,xe,ie] = ode45(@SWFF_RegAdult_NormalSim, tSpan, ics, opts);
% te holds times of sleep onset and wake onset (alternating)
% ie holds event number (1 for sleep and 2 for wake) 
% xe holds value of each of the 7 variables at the event times
% x is vector output of all variables
% t is time

% compute sleep and wake times and durations
sleeponset_og=te(ie==1);
wakeonset_og=te(ie==2);

nbouts=min(length(wakeonset_og),length(sleeponset_og));
sleeponset=sleeponset_og(1:nbouts);
wakeonset=wakeonset_og(1:nbouts);

% calculation duration of sleep and wake bouts
if length(sleeponset_og) > length(wakeonset_og)
    sleepdurs = abs((sleeponset-wakeonset));
    wakedurs = [sleeponset(1); abs((wakeonset_og(1:end)-sleeponset_og(2:end)))];
    extra_S = tSpan(end) - (sleeponset_og(end));
    sleepdurs = [sleepdurs; extra_S];
else
    sleepdurs = abs((sleeponset-wakeonset));
    wakedurs = [sleeponset(1); abs((wakeonset_og(1:end-1)-sleeponset_og(2:end)))];
    extra_W = tSpan(end) - (wakeonset_og(end));
    wakedurs = [wakedurs; extra_W];
end

% for plotting purposes
SleepBlocks = -2*ones(size(t));
for i = 1:length(sleepdurs)
    % put 1 in the entry of the vector when sleep occurs
    start = sleeponset(i); % this is in minutes
    endtime = sleeponset(i) + sleepdurs(i); % this is also in minutes
    SleepBlocks(round(start/dt):round(endtime/dt)) = 0;
end

per=24*60;% minutes -- 24 hours
t2=12*60; % Light duration duration 
Iamp=600; % Lux value -- 
I_vec = NaN(size(t));
for i=1:length(t)
    I_vec(i) = Iamp*heaviside(-1*((t(i)-floor(t(i)/per)*per)-t2));
end
LightBlocks = I_vec == 600;



%% Plot the output of the sleep model
% for x labels in 4 hour increments
% xticknumbers = [0:4:(d*24)]*60;
% xticklabelNames = repmat([8 12 16 20 0 4],1,d+1);

% if want just days instead
xticknumbers = 0:1440:d*1440;
xticklabelNames = string(1:d);


figure 
subplot(3,1,1)
plot(t,x(:,1), 'color',[0.00,0.37,0.74],'LineWidth',3)
hold on
plot(t,x(:,2), 'color',[0.64,0.08,0.18],'LineWidth',3)
% plot light blocks as light yellow
plot(t,15*LightBlocks-5,'*','color',[1.0 1.0 0.6])
% plot sleep blocks
plot(t,SleepBlocks-1,'*','color',[0.6 0.6 0.6])
%rectangle('Position',[sleeponset(1) 0 sleepdurs(1) 25],'facecolor',[0.7 0.7 0.7 0.3],'edgecolor',[0.7 0.7 0.7 0.3])
title('Firing Rate of Sleep Populations')
legend('fWake', 'fSleep')
xticks(xticknumbers)
xticklabels(string(xticklabelNames))
xlabel('Clock time')
ylabel('Firing rates (Hz)')
ylim([-1 10])
xlim([0 d*24*60])
set(gca,'fontsize',18)
subplot(3,1,2)
plot(t,x(:,4),'color',[0.49,0.18,0.56],'LineWidth',3)
hold on
% plot light blocks as light yellow
plot(t,310*LightBlocks-8,'*','color',[1.0 1.0 0.6])
% plot sleep blocks
plot(t,SleepBlocks,'*','color',[0.6 0.6 0.6])
ylim([-1 303])
xlim([0 d*24*60])
xticks(xticknumbers)
xticklabels(string(xticklabelNames))
title('Sleep Homeostatic Drive')
xlabel('Clock time')
ylabel('Slow Wave Activity')
set(gca,'fontsize',18)
subplot(3,1,3)
plot(t,x(:,5),'color',[0.13,0.67,0.64],'LineWidth',3)
hold on
% plot light blocks as light yellow
plot(t,6.2*LightBlocks-5,'*','color',[1.0 1.0 0.6])
% plot sleep blocks
plot(t,SleepBlocks-1.2,'*','color',[0.6 0.6 0.6])
ylim([-1.3 1.2])
xlim([0 d*24*60])
title('Sleep Circadian Component')
xticks(xticknumbers)
xticklabels(string(xticklabelNames))
xlabel('Clock time')
ylabel('Circadian Activity')
set(gca,'fontsize',18)

%% run pain homeostat 

%Plot the circadian transformation from sleep to pain
C_pain = circadian_transformation(x(:,5));
P = (2*pi)/(24*60);

% initial conditions for pain homeostat
swu = -0.45; %initial pain level at onset of wake (from homeostatic)
sso = 0; %initial pain level at onset of sleep
% run pain homeostat
H_pain = pain_homeostat(wakedurs,sleepdurs,wakeonset,sleeponset,t,swu,sso);

%combined pain is sum of circadian and homeostatic drive
P1 = H_pain + C_pain;


% plotting the three components of the model
figure
subplot(3,1,1)
plot(t, P1,'color',[0.00,0.37,0.74],'LineWidth',3)
hold on 
% plot light blocks as light yellow
plot(t,5.9*LightBlocks-5,'*','color',[1.0 1.0 0.6])
% plot sleep blocks
plot(t,SleepBlocks-0.65,'*','color',[0.6 0.6 0.6])
xlabel('Clock Time')
ylabel('Pain Intensity')
title('Total Pain Intensity')
ylim([-0.7 0.9])
xlim([0 d*24*60])
yticks([-.65 -.3 0 .3 .6 .9])
xticks(xticknumbers)
xticklabels(string(xticklabelNames))
set(gca,'fontsize',18)
subplot(3,1,2)
plot(t, H_pain,'color',[0.49,0.18,0.56],'LineWidth',3)
hold on 
% plot light blocks as light yellow
plot(t,5.9*LightBlocks-5,'*','color',[1.0 1.0 0.6])
% plot sleep blocks
plot(t,SleepBlocks-0.65,'*','color',[0.6 0.6 0.6])
xlabel('Clock Time')
ylabel('Pain Intensity')
title('Pain Homeostatic Drive')
ylim([-0.65 0.9])
xlim([0 d*24*60])
yticks([-.6 -.3 0 .3 .6 .9])
xticks(xticknumbers)
xticklabels(string(xticklabelNames))
set(gca,'fontsize',18)
subplot(3,1,3)
plot(t, C_pain,'color',[0.13,0.67,0.64],'LineWidth',3)
hold on 
% plot light blocks as light yellow
plot(t,5.9*LightBlocks-5,'*','color',[1.0 1.0 0.6])
% plot sleep blocks
plot(t,SleepBlocks-0.65,'*','color',[0.6 0.6 0.6])
xlabel('Clock Time')
ylabel('Pain Intensity')
title('Pain Circadian Component')
ylim([-0.65 0.9])
xlim([0 d*24*60])
yticks([-.6 -.3 0 .3 .65 .9])
xticks(xticknumbers)
xticklabels(string(xticklabelNames))
set(gca,'fontsize',18)


function [value, stopflag, direction] = events(~,x)
% events function to save bout onset times
% sleep onset = 1
value(1) = x(1) - 4.0;    %fw decreasing through 4
stopflag(1) = 0;
direction(1) = -1;

% wake onset = 2
value(2) = x(1) - 4.0;    %fw increasing through 4
stopflag(2) = 0;
direction(2) = 1;

end
