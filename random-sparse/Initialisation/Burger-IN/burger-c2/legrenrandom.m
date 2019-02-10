%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       S-SL method for sparse indentification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================================================================
% Target equation: the burger's equation
%       u_t=0.1u_xx+uu_x, 
%
%
% Tunable parameters:
%     lambda:corresopnding the regularization parameter in sparse 
%          optimization algorithm 
%     maxmum_l: The times for the repeated experiment
%     max_num:The running times for each lambda
%     K: Number of burst
%
%
% Subfunctions:
%     dynamics(Define the system dynamics)
%     tac_reconstruction(A iterative re-weighted l1-minimisation algorithm)
%     make_fft(Approximate the deirvatives from data using the Fourier
%            spectal method)
%     leg2mon (Legendre to Monomial Transform)
%     dictionarycombi(Construct the Legrendre dictionary matrix)
%  
%
% Reference: Xiuting Li, Liang Li et al. Sparse Learning of Partial
%          Differential Equations with Structured Dictionary Matrix.
%
%
% Authors: @Xiuting Li
% Date: Jun 22, 2018
% =================================================================

clc; clear all; close all; warning off

% =================================================================
% Steps 1: Generate the data matrix U.
% =================================================================

max_num=1;
lambda=0.6;
for repeat=1:max_num
    lambda=1*lambda;
    K=5;
    for i=1:K
        L=16; n=30;
        x2=linspace(-L/2,L/2,n+1); x=x2(1:n); %space
        k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
        t=[0:5*10e-8/2:5*10e-8]; % time
        dt=t(2)-t(1);
        dx=x(2)-x(1);
        u00=[sech(x);2*sech(x);exp(-(x+2).^2);0.2*exp(-(x+2).^2);3*sech(x);sin(pi*x)];
        for r=1:n
            u(1,r)=u00(i,r);%+0.01*(2*rand(1,1)-1);%
        end
        %u=u00(i,:);%+0.01*(2*rand(1,n)-1); % initial conditions 
        ut=fft(u);
        [t,utsol]=ode45('viscous',t,ut,[],k,x);
        for j=1:length(t)
            usol(j,:)=ifft(utsol(j,:)); % bring back to space
        end
    figure(1)
    mesh(abs(usol));
    data{i} = real(usol([1,end],:)); % record the initial time and final time
    end


%% =============================================================
% Steps 2: Estimating the derivatives.
% =================================================================

    for i=1:K
        U=data{i};
        ut = ( U(2,:) - U(1,:) )/(5*10e-8);  % time derivative
        ut=ut(1,:);   %keep the data at t_0
        Ut{i}=ut.';
    end

x_th = 2;                   % order of derivatives (x_th=2,u_xx)
parameter.deg = 10;         %degree of polynomial to use 
parameter.x_num_to_fit = 5; %number of points to cut for spatial derivatives
parameter.t_num_to_fit = 0; %number of points to cut for time derivative


    for i=1:K
        u0=data{i};  
%         [uleft]  = make_input_poly(u0(1,:).',x,x_th,parameter);
%         Theta{i}=uleft.derivative;
        [uxxxx] = make_fft(u0(1,:).',4,L);
        [uxxx] = make_fft(u0(1,:).',3,L);
        [uxx] = make_fft(u0(1,:).',2,L);
        [ux] = make_fft(u0(1,:).',1,L);
        uleft.derivative=[u0(1,:).',ux,uxx,uxxx,uxxxx];
        %uleft.derivative=derivative(x,i);
        Thetafft{i}=uleft.derivative;
    end


%% =============================================================
% Steps 3: Construct legrendre dictionary matrix.
% =================================================================

% Keep the same length with Theta (due to polynomial interplant)
    for i=1:K
        Utcut{i}=Ut{i}(1:n,1);
        Theta{i}=Thetafft{i}(1:n,:);
    end


% Select the  smaple point
    for i=1:K
        Ycut{i}=Utcut{i}(1:n,:);
        Tcut{i}=Theta{i}(1:n,:);
    end


% Stack
    for i=1:K
        dall((i-1)*(n-1+1)+1:i*(n-1+1),:)=Tcut{i};
        y0((i-1)*(n-1+1)+1:i*(n-1+1),1)=Ycut{i};
    end 
 

%  Scale the data to be valued in [-1, 1].
%  a = 2/( max(dall(:)) - min(dall(:)) );
%  b = -2*min(dall(:))/( max(dall(:)) - min(dall(:)) )-1;
%  dall_1 = a*dall + b;
a=1;b=0;
dall_1 = a*dall + b;
 
 % Compute the true coefficient for 2-order Legrendre basis
%ctrue=[-0.1*b/a-b^2/(a^2), b/(a^2), b/(a^2), 0.1/(a), 0, 0, 0, -1/(a^2), 0, 0];

 % Compute the true coefficient for 4-order derivative
 alpha=100;
ctrue=[-alpha*b/a-b^2/(a^2), b/(a^2), b/(a^2), alpha/(a), 0, 0, 0, 0, 0, 0, 0, -1/(a^2), 0, 0,...
    0,0,0,0,0,0,0];

% Compute the true coefficient for 3-order Legrendre basis
% ctrue=[-0.1*b/a-b^2/(a^2), b/(a^2), b/(a^2), 0.1/(a), 0, 0, 0, -1/(a^2), 0, 0,...
%     0,0,0,0,0,0,0,0,0,0];
    
 
% Test two sides 
right1=b/(a^2)*dall_1(:,1)+b/(a^2)*dall_1(:,2)+0.1/a*dall_1(:,3)...
         -1/(a^2)*dall_1(:,1).*dall_1(:,2)-0.1*b/a-b^2/(a^2);
figure(i+1)
   plot(y0,'LineWidth',2,'Color',[1 0 0])
   hold on
   plot(right1,'MarkerSize',1.5,'Marker','o','LineWidth',1.5,'LineStyle','-.',...
    'Color',[0 0 1])
 

% Construct the Legrendre dictionary matrix
   Tall_1=dictionarycombi(dall_1,2);
  
% Compute the norm of each column.
      Aleg_cnorm = sqrt(sum(Tall_1.^2,1)); 
      Aleg_cnorm = Aleg_cnorm(:);
     % Normalize columns of Aleg.
      Aleg1 = normc(Tall_1);
  
   
%% ========================================================
% Steps 4: Use SBL to identify the coefficients.
% =================================================================
iter=5;
w_estimate = tac_reconstruction( y0, Aleg1,lambda, iter);
for i=1:size(w_estimate,2)
w_estimate(:,i) = w_estimate(:,i)./Aleg_cnorm; 
end

for i=size(w_estimate,2)-1:size(w_estimate,2)
    Ind1=[2, 3, 4,5,6];
    Ind20=[7,8,9,10,11];
    Ind11=[12,13,14,15,16,17,18,19];
    Ind300=[11,12,13];
    Ind210=[14,15,16];
    Ind120=[17,18,19];
    Ind111=[20];
    cM(:,i) = leg2mon(w_estimate(:,i),2,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111);
       disp('identified coefficient');
       disp(cM(:,i)');
       disp('true coefficient');
       disp(ctrue);
       RMS=norm(ctrue(1,:)-cM(:,i)',2)/norm(ctrue(1,:),2);
       disp(RMS);
% err = abs([(ctrue(1,1) -cM(1,i))/ctrue(1,1),  (ctrue(1,2) - cM(2,i))/ctrue(1,2),...
%     (ctrue(1,3) -cM(3,i))/ctrue(1,3),...
%     (ctrue(1,4) -cM(4,i))/ctrue(1,4), (ctrue(1,8) -cM(8,i))/ctrue(1,8)]);
% fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);
end
      cM_all(:,repeat)=cM(:,iter);
      clear cM
end


%% ========================================================
% Steps 5: Stability selection.
% =================================================================
 
[m,n]=size(Tall_1);
index=zeros(n,1);

for i=1:n
    index(i)=size(find(cM_all(i,:)~=0),2);
end

pro=index/max_num;
get_term=index;
threshold=0.7;
for i=1:n
    if pro(i)<threshold
        get_term(i)=0;
    end
end
