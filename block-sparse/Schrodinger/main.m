%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       S-SL method for sparse indentification       %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% Target equation: Schrodinger equation
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
addpath ../Functions
warning off
%% Generate data

% load Time & Space.mat
data = load('harmonic_osc.mat');
u = data.usol;
dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t;
x = data.x;

%% Estimate derivatives

%Method: Pade
t_th = 1; %order of time derivative
x_th = 3; %maximum order of spatial derivative
      
[derivative_r] = make_input_pade(real(u),dx,x_th);  
[derivative_i] = make_input_pade(imag(u),dx,x_th);  

y0_r = make_y_pade(real(u),dt,t_th);
y0_i = make_y_pade(imag(u),dt,t_th);
y0 = y0_r + y0_i*sqrt(-1);

input = derivative_r.derivative+derivative_i.derivative*sqrt(-1);
input_name = derivative_r.name;

%% Builds a dictionary matrix

derdata = input(:,2:end);
dername = input_name(2:end);

data = [input(:,1) reshape(abs(u),size(u,1)*size(u,2),1)];
dataname = {'u','|u|'};

polyorder = 2; %maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(data,derdata,dataname,dername,polyorder);

coeff = zeros(48,290);
yall = zeros(401,300);
A = zeros(80200,4800);

y = zeros(80200,1);

for i = 6:105

    dic = zeros(401,24);
    
    for j = 1:401
        dic(j,:) = theta0((j-1)*300+i,:);
        yall(j,i) = y0((j-1)*300+i);
    end
    
    dic = [real(dic) -imag(dic);imag(dic) real(dic)]; 
    yy = [real(yall(:,i));imag(yall(:,i))];
    
    A((i-6)*802+1:(i-5)*802,(i-6)*48+1:(i-5)*48) = dic;
    y((i-6)*802+1:(i-5)*802) = yy;
        
end

D = sparse(80200,4800);

for i = 1:48
    for j = 1:100
        D(:,j+(i-1)*100) = A(:,(j-1)*48+i);
    end
end

%% Identification

w_ref = zeros(4800,1);
w_ref(3101:3200) = 0.5;
w_ref(2501:2600) = -0.5*x(6:105)'.^2;

partition = 100*ones(48,1);
delta = 1e-4;
lambda = 0.002;
Iter = 10;

w = re_group_lasso(D,y,partition,lambda,delta,Iter);

x_reg = x(6:105)';

theta = [];
for i = 1:7
    theta = [theta x_reg.^(i-1)];
end

theta = [theta sin(x_reg) cos(x_reg) sin(x_reg).^2 cos(x_reg).^2 sin(x_reg).*cos(x_reg)];

lambda2 = 0.005;
lambda3 = 0.00005;

MAXITER = 5;

w1 = tac_reconstruction(w(2501:2600),theta,lambda2,MAXITER);
w1 = w1(:,5);
w1(abs(w1)<1e-6) = 0;

w2 = tac_reconstruction(w(3101:3200),theta,lambda3,MAXITER);
w2 = w2(:,5);
w2(abs(w2)<1e-6) = 0;

%% print result

%RMS Error
fprintf('w1 = %f, RMS error = %f\n',w1(3),norm(w_ref(2501:2600)-w(2501:2600))/norm(w_ref(2501:2600),2));
fprintf('w2 = %f, RMS error = %f\n',w2(1),norm(w_ref(3101:3200)-w(3101:3200))/norm(w_ref(3101:3200),2));
