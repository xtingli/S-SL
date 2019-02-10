%% 
% draw the figures for the identification results
% Target equation: Spatially Advection Diffusion equation: (noise case) 
%       u_t =  (-1.5 + cos(2*pi*x/5))*u_x
%             -2*pi/5*sin(2*pi*x/5)*u
%             +0.1u_xx

%% noisy data

load('SAD_figure.mat'); % data generated from main_noisy.m
% a = sin(x(191:240)).^2;
% b = -x(191:240).^2-x(191:240).^4;
% c = 0.1*ones(1,50);

% u
figure(1)
plot(x(68:167),w_ref(101:200),'LineWidth',3);
hold on 
plot(x(68:167),w(101:200),'LineWidth',3);

% u_x
figure(2)
plot(x(68:167),w_ref(401:500),'LineWidth',3);
hold on 
plot(x(68:167),w(401:500),'LineWidth',3);
axis([-2.5 1.5 -3.0 0]);

% u_xx
figure(3)
plot(x(68:167),w_ref(501:600),'LineWidth',3);
hold on 
plot(x(68:167),w(501:600),'LineWidth',3);
axis([-2.5 1.5 0 0.2]);

