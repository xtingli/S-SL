function [w] = re_group_lasso2(D,y,partition,lambda,delta,MAX_ITER_Reweight)
% re-weighted group lasso
% assume that S = lambda*I

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D: total dictionary matrix based on all C sets of dat
% y: model output vector based on all C sets of data
% partition: the partition of data
% theta:  re-weighted parameter
% lambda: tunning parameter
% delta: threshold for filtering coefficients
%
% Reference: Xiuting Li, Liang Li et al. Sparse Learning of Partial
%          Differential Equations with Structured Dictionary Matrix. 
%
%
% Authors: Xiaoquan Tang, Xiuting Li, Liang Li 
% Date: Feb, 10, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QUIET = 1;
[~, n] = size(D);
ep = 1e-12;
% check that sum(p) = total number of elements in x
if (sum(partition) ~= n)
    error('invalid partition');
end

% the number of groups
Ng = length(partition) ;
% cumulative partition
cum_part = cumsum(partition);
block_idx = cell(Ng);
for k=1:Ng
    block_idx{k} = (cum_part(k)-partition(k)+1:cum_part(k))' ;
end
Dblock = n/Ng;
theta = ones(Ng,1);
alpha = ones(Ng,MAX_ITER_Reweight);
alpha(:,1) = theta/Dblock;
theta_all = ones(Ng,MAX_ITER_Reweight);
gamma = zeros(Ng,MAX_ITER_Reweight);
w_total = [];

for iter=1:1:MAX_ITER_Reweight
    
    if ~QUIET
        fprintf('Iteration %d / total %d iterations.\n', iter, MAX_ITER_Reweight);
    end
    
    % = = = = = = =  = = = = = Solver CVX = = = = = =  = = = = = = = = = = = =
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cvx_begin quiet
    cvx_solver sedumi
    %     cvx_solver SDPT3
    variable v(Dblock,Ng)
    minimize    (lambda*sum(theta.*norms( v, 2 ,1 )')+ 1/2*sum_square(D* vec(v)-y) )
    %                      subject to
    
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w = vec(v);
    %% prune w
%     w(abs(w)./norm(w)<=delta) = 0;
    for i = 1:Ng
        if norm(w((i-1)*partition(1)+1:i*partition(1)))/norm(w) < delta
            w((i-1)*partition+1:i*partition) = 0;
        end
    end
    w_total = [w_total, w];
    
    %% gamma update
    
    for k=1:Ng
        gamma(k,iter)=norm(w(block_idx{k}))/sqrt(Dblock*alpha(k));
    end
    gammas = kron(gamma(:,iter),ones(Dblock,1));
    Gamma = diag(gammas);
    
    %% alpha update    
    Gamma = Gamma + ep*eye(n);
    alpha_new = diag( inv(inv(Gamma)+D'*D) )./(-(gammas+ep).^2) + 1./(gammas+ep);
    
    % theta update
    for k = 1:Ng
        sel = block_idx{k};
        alpha(k,iter) = sum(alpha_new(sel))/length(sel);
    end
    theta = sqrt(abs(Dblock*alpha(:,iter)));
    theta_all(:,iter+1) = theta ;
end

theta_regression = [];
nonzero_index = [];
for i = 1:Ng
    if norm(w((i-1)*partition(1)+1:i*partition(1))) ~= 0
        theta_regression = [theta_regression D(:,(i-1)*partition(1)+1:i*partition(1))];
        nonzero_index = [nonzero_index,i];
    end
end
% 
w_regress = theta_regression\y;
for i = 1:length(nonzero_index)
    w((nonzero_index(i)-1)*partition(1)+1:nonzero_index(i)*partition(1)) = w_regress((i-1)*partition(1)+1:i*partition(1));
end

end