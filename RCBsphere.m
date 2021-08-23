% [as_rE lambda]   =  RCBsphere(R,as_n,epsilon)
%
% INPUT ARGUMENTS
% 'R'   the exact or sample covariance matrix
% 'as_n'    is the nominal steering vector  
% 'epsilon'  is boundary for error mismatch vector, i.e. as_real  =  as_nominal + as_mismatch
%         ||a_mismatch||^2< = epsilon   (for ||as_n||  =  n, where n is number of sensors in array,

%
% OUTPUT ARGUMENT
% 'as_rE ' is the estimate of real steering vector of SOI by RCBsphere method
% 'lambda' is the lagrange multiplier
%
% Robust Adaptive Beamforming Using RCBsphere method 
% copyright Jin Wei; June 28, 2011
% email: jinzhengwei_521@163.com

function  [as_rE lambda]  =  RCBsphere(R,as_n,epsilon)
    [U,D,V] = svd(R);   % SVD composition of R
    gama = diag(D);
    z = U'*as_n;
    %-----------用牛顿迭代法求解lambda----------------------------------
    lambda = 0;
        while abs(norm(abs(z)./(1+lambda*gama))^2-epsilon)>1e-5
            lambda = lambda+(norm(abs(z)./(1+lambda*gama))^2-epsilon)/ ...
                sum(2*abs(z).^2.*gama.*(1+lambda*gama)./(1+lambda*gama).^4);
        end
        as_rE  = as_n-U*diag(1./(1+lambda*gama))*U'*as_n;  %the estimate of real steering vector of SOI by RCBsphere method
