%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       S-SL method for sparse indentification       %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% Target equation: Kuramoto-Sivashinsky equation: (noise case)
%       u_t =  (1 + sin(2*pi*x/20)/4) * uu_x
%             +(-1 + exp(-(x-2)^2/5)/4) * u_xx
%             +(-1 - exp(-(x+2)^2/5)/4) * u_xxxx
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
data = load('KS.mat');
u = data.u;

dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t';
x = data.x;

uclean = u;

rng('default');
rng(0);
u = u + 0.0001*std(reshape(u,1,size(u,1)*size(u,2)))*randn(size(u));

%% Estimate derivatives

%Method: polynomial interpolation
t_th = 1; %order of time derivative
x_th = 4; %maximum order of spatial derivative
      
parameter.deg = 9;           % degree of polynomial to use 
parameter.x_num_to_fit = 20; % number of points to use in polynomial interpolation for spatial derivatives
parameter.t_num_to_fit = 5; % number of points to use in polynomial interpolation for time derivative

[derivative]= make_input_poly(u,x,x_th,parameter);  
parameter.deg = 6;
y0 =  make_y_poly(u,t,t_th,parameter);


input = derivative.derivative;
input_name = derivative.name;

%% Builds a dictionary matrix

polyorder = 3; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(input(:,1),input(:,2:end),input_name(1),input_name(2:end),polyorder);

dic = zeros(512-2*20,512-2*5,20);
for i = 1:20
    dic(:,:,i) = reshape(theta0(:,i),512-2*20,512-2*5);
end

ut = reshape(y0,512-2*20,512-2*5);

x_start_ind = 101-20;
x_end_ind = 400-20;
t_start_ind = 1;
t_end_ind = 250;

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
    norms(i) = max(max(abs(D(:,(i-1)*300+1:i*300))));
    T(:,(i-1)*300+1:i*300) = D(:,(i-1)*300+1:i*300)/norms(i);
end

D = sparse(D);
T = sparse(T);

uu_x_true = 1 + 0.25*sin(x*2*pi/20);
u_xx_true = -1 + 0.25*exp(-(x-2).^2/5);
u_4x_true = -1 - 0.25*exp(-(x+2).^2/5);

%% Identification

w_ref = zeros((x_end_ind-x_start_ind+1)*20,1);
w_ref(8*(x_end_ind-x_start_ind+1)+1:9*(x_end_ind-x_start_ind+1)) = uu_x_true(x_start_ind+20:x_end_ind+20)';
w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)) = u_xx_true(x_start_ind+20:x_end_ind+20)';
w_ref(7*(x_end_ind-x_start_ind+1)+1:8*(x_end_ind-x_start_ind+1)) = u_4x_true(x_start_ind+20:x_end_ind+20)';

partition = (x_end_ind-x_start_ind+1)*ones(20,1);
delta = 0.15;
lambda = 0.01;
Iter = 20;

w = re_group_lasso(T, y,partition,lambda,delta,Iter);
for i = 1:20
    w((i-1)*300+1:i*300) = w((i-1)*300+1:i*300)/norms(i);
end

%% print result

fprintf('RMS error (u*u_x) = %f\n',norm(w_ref(8*(x_end_ind-x_start_ind+1)+1:9*(x_end_ind-x_start_ind+1))-w(8*(x_end_ind-x_start_ind+1)+1:9*(x_end_ind-x_start_ind+1)))...
    /norm(w_ref(8*(x_end_ind-x_start_ind+1)+1:9*(x_end_ind-x_start_ind+1)),2));
fprintf('RMS error  (u_xx) = %f\n',norm(w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1))-w(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)))...
    /norm(w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)),2));
fprintf('RMS error (u_xxxx) = %f\n',norm(w_ref(7*(x_end_ind-x_start_ind+1)+1:8*(x_end_ind-x_start_ind+1))-w(7*(x_end_ind-x_start_ind+1)+1:8*(x_end_ind-x_start_ind+1)))...
    /norm(w_ref(7*(x_end_ind-x_start_ind+1)+1:8*(x_end_ind-x_start_ind+1)),2));
