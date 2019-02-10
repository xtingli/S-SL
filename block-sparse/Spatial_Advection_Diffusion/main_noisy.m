%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       S-SL method for sparse indentification       %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% Target equation: Spatially Advection Diffusion equation: (noise case) 
%       u_t =  (-1.5 + cos(2*pi*x/5))*u_x
%             -2*pi/5*sin(2*pi*x/5)*u
%             +0.1u_xx
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
% Data generated from:
% https://github.com/snagcliffs/parametric-discovery
% Rudy, Samuel, et al. "Data-driven identification of parametric partial 
% differential equations." arXiv preprint arXiv:1806.00732 (2018)
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
data = load('Spatially_Advection_Diffusion.mat');
u = data.u;

dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t';
x = data.x;

uclean = u;

rng('default');
rng(0);
u = u + 0.01*std(reshape(u,1,size(u,1)*size(u,2)))*randn(size(u));

%% Estimate derivatives

%Method: polynomial interpolation
t_th = 1; %order of time derivative
x_th = 4; %maximum order of spatial derivative
      
parameter.deg = 9;           % degree of polynomial to use 
parameter.x_num_to_fit = 22; % number of points to use in polynomial interpolation for spatial derivatives
parameter.t_num_to_fit = 14; % number of points to use in polynomial interpolation for time derivative

[derivative]= make_input_poly(u,x,x_th,parameter);  
parameter.deg = 4;
y0 =  make_y_poly(u,t,t_th,parameter);

input = derivative.derivative;
input_name = derivative.name;

%% Builds a dictionary matrix

polyorder = 3; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(input(:,1),input(:,2:end),input_name(1),input_name(2:end),polyorder);

dic = zeros(256-2*22,256-2*14,20);
for i = 1:20
    dic(:,:,i) = reshape(theta0(:,i),256-2*22,256-2*14);
end

ut = reshape(y0,256-2*22,256-2*14);

x_start_ind = 68-22;
x_end_ind = 167-22;
t_start_ind = 1;
t_end_ind = 228;

D = zeros((x_end_ind-x_start_ind+1)*(t_end_ind-t_start_ind+1),(x_end_ind-x_start_ind+1)*20);
y = zeros((x_end_ind-x_start_ind+1)*(t_end_ind-t_start_ind+1),1);

for i = 1:20
    for j = 1:(x_end_ind-x_start_ind+1)
        y((j-1)*(t_end_ind-t_start_ind+1)+1:j*(t_end_ind-t_start_ind+1)) = ut(x_start_ind+j-1,t_start_ind:t_end_ind)';
        D((j-1)*(t_end_ind-t_start_ind+1)+1:j*(t_end_ind-t_start_ind+1),(i-1)*(x_end_ind-x_start_ind+1)+j) = dic(x_start_ind+j-1,t_start_ind:t_end_ind,i)';
    end
end

% normalization
norms = zeros(1,20);
T = zeros(size(D));

for i = 1:20
    norms(i) = max(max(abs(D(:,(i-1)*100+1:i*100))));
    T(:,(i-1)*100+1:i*100) = D(:,(i-1)*100+1:i*100)/norms(i);
end

u_x_true = -1.5 + 1.0*cos(2*x*pi/5);
u_true = -2*pi/5*sin(2*x*pi/5);

w_ref = zeros((x_end_ind-x_start_ind+1)*20,1);
w_ref(4*(x_end_ind-x_start_ind+1)+1:5*(x_end_ind-x_start_ind+1)) = u_x_true(x_start_ind+22:x_end_ind+22)';
w_ref(1*(x_end_ind-x_start_ind+1)+1:2*(x_end_ind-x_start_ind+1)) = u_true(x_start_ind+22:x_end_ind+22)';
w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)) = 0.1;

%% Identification

partition = (x_end_ind-x_start_ind+1)*ones(20,1);
delta = 0.15;
lambda = 0.021;
Iter = 20;

w = re_group_lasso(T, y,partition,lambda,delta,Iter);

for i = 1:20
    w((i-1)*100+1:i*100) = w((i-1)*100+1:i*100)/norms(i);
end

%% print result

fprintf('RMS error (u_x) = %f\n',norm(w_ref(4*(x_end_ind-x_start_ind+1)+1:5*(x_end_ind-x_start_ind+1))-w(4*(x_end_ind-x_start_ind+1)+1:5*(x_end_ind-x_start_ind+1)))...
    /norm(w_ref(4*(x_end_ind-x_start_ind+1)+1:5*(x_end_ind-x_start_ind+1)),2));
fprintf('RMS error  (u) = %f\n',norm(w_ref(1*(x_end_ind-x_start_ind+1)+1:2*(x_end_ind-x_start_ind+1))-w(1*(x_end_ind-x_start_ind+1)+1:2*(x_end_ind-x_start_ind+1)))...
    /norm(w_ref(1*(x_end_ind-x_start_ind+1)+1:2*(x_end_ind-x_start_ind+1)),2));
fprintf('RMS error (u_xx) = %f\n',norm(w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1))-w(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)))...
    /norm(w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)),2));
