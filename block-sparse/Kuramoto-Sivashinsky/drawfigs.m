%% 
% draw the figures for the identification results
% Target equation: Kuramoto-Sivashinsky equation: (noise case)
%       u_t =  (1 + sin(2*pi*x/20)/4) * uu_x
%             +(-1 + exp(-(x-2)^2/5)/4) * u_xx
%             +(-1 - exp(-(x+2)^2/5)/4) * u_xxxx
%% noisy data

load('KS_figure.mat'); % data generated from main_noisy.m

% u*u_x
figure(1)
plot(x(101:400),w_ref(2401:2700),'LineWidth',3);
hold on 
plot(x(101:400),w(2401:2700),'LineWidth',3);
axis([-15 15 0 1.5]);

% u_xx
figure(2)
plot(x(101:400),w_ref(1501:1800),'LineWidth',3);
hold on 
plot(x(101:400),w(1501:1800),'LineWidth',3);
axis([-15 15 -1.5 0]);

% u_xxxx
figure(3)
plot(x(101:400),w_ref(2101:2400),'LineWidth',3);
hold on 
plot(x(101:400),w(2101:2400),'LineWidth',3);
axis([-15 15 -1.5 0]);

