%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       S-SL method for sparse indentification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================================================================
% Target equation: the burger's equation
%       u_t=0.1u_xx+uu_x, u(x,0)=exp(-(x+2)^2)
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

lambda=0.1;
maxmum_l=10;
for repeat0=1:maxmum_l
 lambda=1.01*lambda;
 max_num=20;
for repeat=1:max_num
    K=6;
    for i=1:K
        L=16; n=60;
        x2=linspace(-L/2,L/2,n+1); x=x2(1:n); %space
        k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
        t=[0:5*10e-8/2:5*10e-8]; % time
        dt=t(2)-t(1);
        dx=x(2)-x(1);
        u00=[sech(x);cos(x);exp(-(x+2).^2);3*sech(x);0.2*exp(-(x+2).^2); sin(pi*x);...
            0.3*exp(-(x+2).^2); 2*sin(pi*x);0.5*sin(pi*x);2*sech(x);sech(x);sech(x)];
        for r=1:n
            u(1,r)=u00(i,r)+0.02*(2*rand(1,1)-1);%
        end
        %u=u00(i,:);%+0.01*(2*rand(1,n)-1); % initial conditions 
        ut=fft(u);
        [t,utsol]=ode45('dynamics',t,ut,[],k,x);
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
        
        [uxx] = make_fft(u0(1,:).',2,L);
        [ux] = make_fft(u0(1,:).',1,L);
        uleft.derivative=[u0(1,:).',ux,uxx];
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
 

 % Scale the data to be valued in [-1, 1].
 a = 2/( max(dall(:)) - min(dall(:)) );
 b = -2*min(dall(:))/( max(dall(:)) - min(dall(:)) )-1;
 dall_1 = a*dall + b;
 
 % Compute the true coefficient for 2-order Legrendre basis
 alpha=0.1;
ctrue=[-alpha*b/a-b^2/(a^2), b/(a^2), b/(a^2), alpha/(a), 0, 0, 0, -1/(a^2), 0, 0];


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

w_estimate = tac_reconstruction( y0, Aleg1,lambda, 8);
cLS=pinv([Aleg1(:,1:4),Aleg1(:,8)])*y0;
cLS=[cLS(1),cLS(2),cLS(3),cLS(4),0,0,0,cLS(5),0,0]';
cLS=cLS./Aleg_cnorm;
for i=1:size(w_estimate,2)
w_estimate(:,i) = w_estimate(:,i)./Aleg_cnorm; 
end 

for i=size(w_estimate,2)-1:size(w_estimate,2)
    Ind1=[2, 3, 4];
    Ind20=[5,6,7];
    Ind11=[8,9,10];
    Ind300=[11,12,13];
    Ind210=[14,15,16];
    Ind120=[17,18,19];
    Ind111=[20];
    cM(:,i) = leg2mon(w_estimate(:,i),2,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111);
    cLSM = leg2mon(cLS,2,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111);
    cLStrue=[cLSM(1)+b*cLSM(2)+b*cLSM(3)+b*cLSM(4)+b.^2*cLSM(8),a*cLSM(2)+a*b*cLSM(8),a*cLSM(3)+a*b*cLSM(8),a*cLSM(4),a*a*cLSM(8)];
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
      cM_all(:,repeat)=cM(:,8);
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
threshold=0.6;
for i=1:n
    if pro(i)<threshold
        get_term(i)=0;
    end
end

get_term_all(:,repeat0)=get_term;
pro_all(:,repeat0)=pro;
clear get_term


end

%% max probability of stable variable
for i=1:n
    index_ss(i,:)=max(get_term_all(i,:));
    index_ss_pro(i,:)=max(get_term_all(i,:))/max_num;
end

%% plot 
figure
for i=1:n
plot(pro_all(i,:));
hold on
end
