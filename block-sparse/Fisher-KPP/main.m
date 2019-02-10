%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       S-SL method for sparse indentification       %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% Target equation: Fisher-KPP equation
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
addpath ../Functions
warning off
%% Generate data

% load Time & Space.mat
data = load('fisher.mat');
u = data.u;
dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t;
x = data.x;


%% Estimate derivatives

%Method: Pade
t_th = 1; %order of time derivative
x_th = 3; %maximum order of spatial derivative
      
[derivative] = make_input_pade( u,dx,x_th);    

input = derivative.derivative;
input_name = derivative.name;

y0 = make_y_pade(u,dt,t_th);

%% Builds a dictionary matrix

polyorder = 3; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(input(:,1),input(:,2:end),input_name(1),input_name(2:end),polyorder);

dic = zeros(512,1000,16);
for i = 1:16
    dic(:,:,i) = reshape(theta0(:,i),512,1000);
end

ut = reshape(y0,512,1000);

x_start_ind = 191;
x_end_ind = 240;
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

[U,S,v] = svd(T);
T = U(:,1:800)'*T;
y = U(:,1:800)'*y;
clear U;


%% Identification

partition = (x_end_ind-x_start_ind+1)*ones(16,1);
delta = 1e-1;
Iter = 20;
lambda = 6e-8;

w = re_group_lasso(T, y,partition,lambda,delta,Iter);
for i = 1:16
    w((i-1)*50+1:i*50) = w((i-1)*50+1:i*50)/norms(i);
end

w_ref = zeros((x_end_ind-x_start_ind+1)*16,1);
w_ref(x_end_ind-x_start_ind+2:2*(x_end_ind-x_start_ind+1)) = sin(x(x_start_ind:x_end_ind)').^2;
w_ref(2*(x_end_ind-x_start_ind+1)+1:3*(x_end_ind-x_start_ind+1)) = -(x(x_start_ind:x_end_ind)'.^2+x(x_start_ind:x_end_ind)'.^4);
w_ref(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1)) = 0.1;



x_reg = x(x_start_ind:x_end_ind)';

y1 = w(x_end_ind-x_start_ind+2:2*(x_end_ind-x_start_ind+1),1);
y2 = w(2*(x_end_ind-x_start_ind+1)+1:3*(x_end_ind-x_start_ind+1),1);
y3 = w(5*(x_end_ind-x_start_ind+1)+1:6*(x_end_ind-x_start_ind+1),1);

theta = [];
for i = 1:9
    theta = [theta x_reg.^(i-1)]; 
end

theta = [theta sin(x_reg) cos(x_reg) sin(x_reg).^2 cos(x_reg).^2 sin(x_reg).*cos(x_reg)];

MAXITER = 10;

lambda1 = 6e-5;
lambda2 = 0.0077;
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

