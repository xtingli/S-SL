%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       S-SL method for sparse indentification       %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% Target equation: Fisher-KPP equation (noise case)
%       u_t = 0.1u_xx + sin(x)^2*u - (x^2 + x^4)*u^2
%
%
% Tunable parameters:
%     lambda: corresopnding the regularization parameter in sparse 
%          optimization algorithm 
%     delta: threshold for filtering coefficients
%     Iter: maximum number of iteration 
%
%
% Reference: Xiuting Li, Liang Li et al. Sparse Learning of Partial
%          Differential Equations with Structured Dictionary Matrix. 
%
%
% Author: Liang Li (lli@alumni.hust.edu.cn)
% Date: Feb, 10, 2019
% =========================================================================

% clc
% clear 
% close all
warning off
addpath ../Functions
%% Generate data

% load Time & Space.mat
data = load('fisher.mat');
u = data.u;
dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t;
x = data.x;

t = t(1:800);
x = x(171:260);

u = u(171:260,1:800);

rng('default');
rng(0);
u = u + 0.00001*std(reshape(u,1,size(u,1)*size(u,2)))*randn(size(u));

%% Estimate derivatives

% Method: polynomial interpolation/Pade
t_th = 1; %order of time derivative
x_th = 3; %maximum order of spatial derivative
      
parameter.deg = 9;           % degree of polynomial to use 
parameter.x_num_to_fit = 10; % number of points to use in polynomial interpolation for spatial derivatives
parameter.t_num_to_fit = 0;  % number of points to use in polynomial interpolation for time derivative

[derivative]= make_input_poly(u,x,x_th,parameter);  

y0 = make_y_pade(u,dt,t_th);
ut = reshape(y0,90,800);
ut = ut(11:end-10,:);

input = derivative.derivative;
input_name = derivative.name;

%% Builds a dictionary matrix

polyorder = 3; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(input(:,1),input(:,2:end),input_name(1),input_name(2:end),polyorder);

dic = zeros(90-2*parameter.x_num_to_fit,800-2*parameter.t_num_to_fit,16);
for i = 1:16
    dic(:,:,i) = reshape(theta0(:,i),90-2*parameter.x_num_to_fit,800-2*parameter.t_num_to_fit);
end

x_start_ind = 21-10;
x_end_ind = 70-10;
t_start_ind = 8;
t_end_ind = 407;

D = zeros((x_end_ind-x_start_ind+1)*(t_end_ind-t_start_ind+1),(x_end_ind-x_start_ind+1)*16);
y = zeros((x_end_ind-x_start_ind+1)*(t_end_ind-t_start_ind+1),1);

for i = 1:16
    for j = 1:(x_end_ind-x_start_ind+1)
        y((j-1)*(t_end_ind-t_start_ind+1)+1:j*(t_end_ind-t_start_ind+1)) = ut(x_start_ind+j-1,t_start_ind:t_end_ind)';
        D((j-1)*(t_end_ind-t_start_ind+1)+1:j*(t_end_ind-t_start_ind+1),(i-1)*(x_end_ind-x_start_ind+1)+j) = dic(x_start_ind+j-1,t_start_ind:t_end_ind,i)';
    end
end

% normalization
norms = zeros(1,16);
T = zeros(size(D));

for i = 1:16
    norms(i) = max(max(abs(D(:,(i-1)*50+1:i*50))));
    T(:,(i-1)*50+1:i*50) = D(:,(i-1)*50+1:i*50)/norms(i);
end

w_ref = zeros((x_end_ind-x_start_ind+1)*16,1);
w_ref(x_end_ind-x_start_ind+2:2*(x_end_ind-x_start_ind+1)) = sin(x(x_start_ind+parameter.x_num_to_fit:x_end_ind+parameter.x_num_to_fit)').^2;
w_ref(2*(x_end_ind-x_start_ind+1)+1:3*(x_end_ind-x_start_ind+1)) = -(x(x_start_ind+parameter.x_num_to_fit:x_end_ind+parameter.x_num_to_fit)'.^2+x(x_start_ind+parameter.x_num_to_fit:x_end_ind+parameter.x_num_to_fit)'.^4);
w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)) = 0.1;

%% Identification

partition = (x_end_ind-x_start_ind+1)*ones(16,1);
lambda = 2e-5;
delta = 0.04;
Iter = 30;

w = re_group_lasso2(T,y,partition,lambda,delta,Iter);
for i = 1:16
    w((i-1)*50+1:i*50) = w((i-1)*50+1:i*50)/norms(i);
end

x_reg = x(x_start_ind+10:x_end_ind+10)';

y1 = w(x_end_ind-x_start_ind+2:2*(x_end_ind-x_start_ind+1),1);
y2 = w(2*(x_end_ind-x_start_ind+1)+1:3*(x_end_ind-x_start_ind+1),1);
y3 = w(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1),1);

theta = [];
for i = 1:9
    theta = [theta x_reg.^(i-1)]; 
end

theta = [theta sin(x_reg) cos(x_reg) sin(x_reg).^2 cos(x_reg).^2 sin(x_reg).*cos(x_reg)];

MAXITER = 10;

lambda1 = 8.9*1e-8;
lambda2 = 9.1*1e-6;
lambda3 = 4e-5;

w1 = tac_reconstruction(y1,theta,lambda1,MAXITER);
w1 = w1(:,end);
w2 = tac_reconstruction(y2,theta,lambda2,MAXITER);
w2 = w2(:,end);
w3 = tac_reconstruction(y3,theta,lambda3,MAXITER);
w3 = w3(:,end);

%% print result

%RMS Error
fprintf('w1 = %f, RMS error = %f\n',w1(12),norm(w_ref(x_end_ind-x_start_ind+2:2*(x_end_ind-x_start_ind+1))-w(x_end_ind-x_start_ind+2:2*(x_end_ind-x_start_ind+1),:))...
    /norm(w_ref(x_end_ind-x_start_ind+2:2*(x_end_ind-x_start_ind+1)),2));
fprintf('w2 = %f, %f, RMS error = %f\n',w2(3),w2(5),norm(w_ref(2*(x_end_ind-x_start_ind+1)+1:3*(x_end_ind-x_start_ind+1))-w(2*(x_end_ind-x_start_ind+1)+1:3*(x_end_ind-x_start_ind+1)))...
    /norm(w_ref(2*(x_end_ind-x_start_ind+1)+1:3*(x_end_ind-x_start_ind+1)),2));
fprintf('w3 = %f, RMS error = %f\n',w3(1),norm(w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1))-w(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)))...
    /norm(w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)),2));


