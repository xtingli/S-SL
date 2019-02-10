%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       S-SL method for sparse indentification       %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% Target equation: Schrodinger equation: (noise case) 
%       u_t = 0.5iu_xx - 0.5ix^2*u
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


clc
clear 
close all
warning off
addpath ../Functions
%% Generate data

% load Time & Space.mat
data = load('harmonic_osc.mat');
u = data.usol;
dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t;
x = data.x;

rng('default');
rng(0);
u = u + 0.01/sqrt(2)*std(reshape(real(u),1,size(u,1)*size(u,2)))*randn(size(u))+0.01/sqrt(2)*std(reshape(imag(u),1,size(u,1)*size(u,2)))*randn(size(u));

%% Estimate derivatives

%Method: polynomial interpolation
t_th = 1; %order of time derivative
x_th = 3; %maximum order of spatial derivative
     
parameter.deg = 6;           %degree of polynomial to use 
parameter.x_num_to_fit = 20; %number of points to use in polynomial interpolation for spatial derivatives
parameter.t_num_to_fit = 10; %number of points to use in polynomial interpolation for time derivative

[derivative_r]= make_input_poly( real(u),x,x_th,parameter);
[derivative_i]= make_input_poly( imag(u),x,x_th,parameter);

y0_r =  make_y_poly( real(u),t,t_th,parameter );     
y0_i =  make_y_poly( imag(u),t,t_th,parameter );  
y0 = y0_r + y0_i*sqrt(-1);   

input = derivative_r.derivative+derivative_i.derivative*sqrt(-1);
input_name = derivative_r.name;

%% Builds a dictionary matrix

x_num = parameter.x_num_to_fit;
t_num = parameter.t_num_to_fit;

derdata = input(:,2:end);
dername = input_name(2:end);

data = [input(:,1) reshape(abs(u(x_num+1:end-x_num,t_num+1:end-t_num)),(size(u,1)-2*x_num)*(size(u,2)-2*t_num),1)];
dataname = {'u','|u|'};

polyorder = 2; %maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(data,derdata,dataname,dername,polyorder);

coeff = zeros(48,290);
yall = zeros(381,260);
A = zeros(38100,2400);

y = zeros(38100,1);

for i = 61:110

    dic = zeros(381,24);
    
    for j = 1:381
        dic(j,:) = theta0((j-1)*260+i,:);
        yall(j,i) = y0((j-1)*260+i);
    end
    
    dic = [real(dic) -imag(dic);imag(dic) real(dic)]; 
    yy = [real(yall(:,i));imag(yall(:,i))];
    
    A((i-61)*762+1:(i-60)*762,(i-61)*48+1:(i-60)*48) = dic;
    y((i-61)*762+1:(i-60)*762) = yy;
        
end

D = sparse(38100,2400);

for i = 1:48
    for j = 1:50
        D(:,j+(i-1)*50) = A(:,(j-1)*48+i);
    end
end

w_ref = zeros(2400,1);
w_ref(1551:1600) = 0.5;
w_ref(1251:1300) = -0.5*x(81:130)'.^2;

%% Identification

partition = 50*ones(48,1);
delta = 1e-4;
lambda = 0.07;
Iter = 10;

w = re_group_lasso(D, y,partition,lambda,delta,Iter);

x_reg = x(81:130)';

theta = [];
for i = 1:7
    theta = [theta x_reg.^(i-1)];
end

theta = [theta sin(x_reg) cos(x_reg) sin(x_reg).^2 cos(x_reg).^2 sin(x_reg).*cos(x_reg)];

lambda2 = 0.12;
lambda3 = 0.000075;

MAXITER = 5;

w1 = tac_reconstruction(w(1251:1300),theta,lambda2,MAXITER);
w1 = w1(:,5);
w1(abs(w1)<1e-6) = 0;

w2 = tac_reconstruction(w(1551:1600),theta,lambda3,MAXITER);
w2 = w2(:,5);
w2(abs(w2)<1e-6) = 0;

%% print result

fprintf('w1 = %f, RMS error = %f\n',w1(3),norm(w_ref(1251:1300)-w(1251:1300))/norm(w_ref(1251:1300),2));
fprintf('w2 = %f, RMS error = %f\n',w2(1),norm(w_ref(1551:1600)-w(1551:1600))/norm(w_ref(1551:1600),2));
