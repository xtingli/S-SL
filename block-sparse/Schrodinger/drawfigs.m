%% 
% draw the figures for the identification results
% Target equation: Schrodinger equation
%       u_t = 0.5iu_xx - 0.5ix^2*u


%% clean data

% load('drawfig.mat'); % data generated from main.m
% b = (-x(6:105).^2/2)';
% a = 0.5*ones(100,1);
% 
% figure(1)
% plot(x(6:105),b,'LineWidth',3);
% hold on 
% plot(x(6:105),w(2501:2600),'LineWidth',3);
% 
% figure(2)
% plot(x(6:105),a,'LineWidth',3);
% hold on 
% plot(x(6:105),w(3101:3200),'LineWidth',3);
% axis([-4.5 -1 0 0.6]);



%% noisy data 

load('drawfig_noisy.mat'); % data generated from main_noisy.m
b = (-x(81:130).^2/2)';
a = 0.5*ones(50,1);

figure(1)
plot(x(81:130),b,'LineWidth',3);
hold on 
plot(x(81:130),w(1251:1300),'LineWidth',3);

figure(2)
plot(x(81:130),a,'LineWidth',3);
hold on 
plot(x(81:130),w(1551:1600),'LineWidth',3);
axis([-2.2 -1 0 0.6]);