%% 
% draw the figures for the identification results
% Target equation: Fisher-KPP equation
%       u_t = 0.1u_xx + sin(x)^2*u - (x^2 + x^4)*u^2

%% clean data
% 
% load('fisherfig.mat'); % data generated from main.m
% a = sin(x(191:240)).^2;
% b = -x(191:240).^2-x(191:240).^4;
% c = 0.1*ones(1,50);
% 
% 
% figure(1)
% plot(x(191:240),a,'LineWidth',3);
% hold on 
% plot(x(191:240),w(51:100),'LineWidth',3);
% axis([-1.6 -0.2 0 1]);
% 
% figure(2)
% plot(x(191:240),b,'LineWidth',3);
% hold on 
% plot(x(191:240),w(101:150),'LineWidth',3);
% 
% figure(3)
% plot(x(191:240),c,'LineWidth',3);
% hold on 
% plot(x(191:240),w(251:300),'LineWidth',3);
% axis([-1.6 -0.2 0 0.12]);

%% noisy data 
load('fisherfig_noisy.mat'); %data generated from main_noisy.m
a = sin(x(191:240)).^2;
b = -x(191:240).^2-x(191:240).^4;
c = 0.1*ones(1,50);


figure(1)
plot(x(191:240),a,'LineWidth',3);
hold on 
plot(x(191:240),w(51:100),'LineWidth',3);
axis([-1.6 -0.2 0 1]);

figure(2)
plot(x(191:240),b,'LineWidth',3);
hold on 
plot(x(191:240),w(101:150),'LineWidth',3);

figure(3)
plot(x(191:240),c,'LineWidth',3);
hold on 
plot(x(191:240),w(251:300),'LineWidth',3);
axis([-1.6 -0.2 0 0.12]);
