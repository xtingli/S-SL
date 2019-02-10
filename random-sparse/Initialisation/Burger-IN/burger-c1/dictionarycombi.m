function phiX = dictionarycombi(U, order)
% Description: Construct the dictionary matrix phiX containing all
% multivariate monomials up to degree two or three
% Input: U = [u(x1,t1) ux(x1,t1) ... uxxx(x1,t1)
%             u(x2,t1) ux(x2,t1) ... uxxx(x2,t1)
%                    ......
%             u(xn,tm) ux(xn,tm) ... uxxx(xn,tm)]
%        order = degree of (monomial) or 'legendre'
% Output: the dictionary matrix phiX of size mn by N, where m= #measurements
% and N = total numnber of the candidate terms
% Author:Xiuting Li.
% Date: Jun 22, 2018

if (order<=2)
m = size(U,1); % number of measurements
n = size(U,2); % number of condidate term
p = order;

phiX = zeros(m,factorial(n+p)/(factorial(n)*factorial(p)));

% 1 - 1 column
U0(:,1) = ones(m,1);   % 1

% X - n columns - terms of the form sqrt(3)*ui
U1(:,1:n)=sqrt(3)*U;   % x1 x2 ... xn


% X^2 - n*(n+1)/2 columns - terms of the form (sqrt(5)/2.0)*(3*(ui)^2-1)
for k = 1:n
    U2(:,k)=(sqrt(5)/2)*(3*U(:,k).*U(:,k)-ones(m,1));
end


% X^3 -  terms of the form 3*(ui)*(uj)
ind3=1;
for k = 1:n-1
    for i=k+1:n
        U3(:,ind3)=3*U(:,k).*U(:,i);
        ind3 = ind3+1;
    end
end
phiX=[U0, U1, U2, U3];
end


if (order<=3)&&(order>=3)
   m = size(U,1); % number of measurements
   n = size(U,2); % dimension of the ODE
   p = order;
   
   phiX = zeros(m,factorial(n+p)/(factorial(n)*factorial(p)));% 1 + n + n*(n+1)/2

% 1 - 1 column
U0(:,1) = ones(m,1);   % 1

% X - n columns - terms of the form sqrt(3)*ui
U1(:,1:n)=sqrt(3)*U;   % x1 x2 ... xn


% X^2 - n*(n+1)/2 columns - terms of the form (sqrt(5)/2.0)*(3*(ui)^2-1)
for k = 1:n
    U2(:,k)=(sqrt(5)/2)*(3*U(:,k).*U(:,k)-ones(m,1));
end


% X^3 -  terms of the form 3*(ui)*(uj)
ind3=1;
for k = 1:n-1
    for i=k+1:n
        U3(:,ind3)=3*U(:,k).*U(:,i);
        ind3 = ind3+1;
    end
end
% terms of the form sqrt(7)/2*(5*(ui)^3-ui)
for k = 1:n
    U4(:,k) = (sqrt(7)/2)*(5*U(:,k).^3-3*U(:,k));
end


% terms of the form sqrt(15)/2*(3*(ui)^2-1)*(uk)
ind5=1;
for k = 1:n-1
    for i=k+1:n
        U5(:,ind5)=(sqrt(15)/2)*(3*U(:,k).^2-ones(m,1)).*U(:,i);
        ind5 = ind5+1;
    end
end

ind6=1;
for k = 1:n-1
    for i = k+1:n
        U6(:,ind6)= (sqrt(15)/2)*(3*U(:,i).^2-ones(m,1)).*U(:,k);
        ind6 = ind6+1;
    end
end

 
% terms of the form sqrt(27)*(ui)*(uj)*(uk)
ind7=1;
for k = 1:n-2
    for i=k+1:n-1
        for j=i+1:n
            U7(:, ind7) = sqrt(27)*U(:,k).*U(:, i).*U(:, j);
            ind7 = ind7+1;
        end
    end
end
phiX=[U0, U1, U2, U3, U4, U5, U6, U7];
end
