function cM = leg2mon(cL,p,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111)

% ============================================================
% Inputs:
%   cL = coefficient wrt the dictionary which consists of
%           Legendre polynomials of degree at most 3
%   p = maximum degree of polynomials in the dictionary
%
%   Ind are the indices for each of the monomial terms:
%
%   Ind1 = indices for terms of the form ui
%   Ind20 = indices for terms of the form (ui)^2
%   Ind11 = indices for terms of the form (ui)*(uj)
%   Ind300 = indices for terms of the form (ui)^3
%   Ind210 = indices for terms of the form (ui)^2*(uj)
%   Ind120 = indices for terms of the form (ui)*(uj)^2
%   Ind111 = indices for terms of the form (ui)*(uj)*(uk)
%
% Output:
%   cM = coefficient wrt the dictionary which consists of
%           monomials of degree at most 3
%
% Recall: Legendre polynomials
%   P0: 1
%   P1: sqrt(3)*ui
%   P2: (3*ui^2-1)*sqrt(5)/2
%       3*ui*uj
%   P3: (5*ui^3-3*ui)*sqrt(7)/2
%       sqrt(15)*(3*ui^2-1)*uj
%       sqrt(27)*ui*uj*uk
%
% Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
% Date: May 9, 2018
% ============================================================


% initialization
cM = zeros(size(cL));

% 1
cM(1) = cL(1)- sum(cL(Ind20))*sqrt(5)/2; %constant

if p==2
    
    % ui
    cM(Ind1) = cL(Ind1) * sqrt(3);
    %(ui)^2
    cM(Ind20) = cL(Ind20) * sqrt(5) * 3/2;
    % (ui)*(uj)
    cM(Ind11) = cL(Ind11) * 3;
    
elseif p==3
    
    % ui
    %   terms from Ind1 and Ind300
    cM(Ind1) = cL(Ind1) * sqrt(3) - cL(Ind300) * sqrt(7) * 3/2;
    n = length(Ind1);
    %   terms from Ind210
    s = 1;
    for ii=2:n
        cM(Ind1(ii:n)) = cM(Ind1(ii:n)) - cL(Ind210(s:s+n-ii)) * sqrt(15)/2;
        s = s + n - ii + 1;
    end
    %   terms from Ind120
    s = 1;
    for ii=1:n-1
        cM(Ind1(ii)) = cM(Ind1(ii)) - sum(cL(Ind120(s:s+n-ii-1))) * sqrt(15)/2;
        s = s + n - ii;
    end
    
    %(ui)^2
    cM(Ind20) = cL(Ind20) * sqrt(5) * 3/2;
    % (ui)*(uj)
    cM(Ind11) = cL(Ind11) * 3;
    
    %(ui)^3
    cM(Ind300) = cL(Ind300) * sqrt(7)*5/2;
    % (ui)^2*(uj)
    cM(Ind210) = cL(Ind210) * sqrt(15)*3/2;
    % (ui)*(uj)^2
    cM(Ind120) = cL(Ind120) * sqrt(15)*3/2;
    % (ui)*(uj)*(uk)
    cM(Ind111) = cL(Ind111) * sqrt(3)*3;
    
end
end