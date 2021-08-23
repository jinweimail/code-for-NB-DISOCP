% Robust Beamforming Using LCSS
% copyright Guo Yunzhou; Aug 25, 2020
% email: mailguoyunzhou@163.com

% according to paper "AMAR A, DORON MA. A linearly constrained minimum variance 
%                        beamformer with a pre-specified suppression level over
%                        a pre-defined broad null sector[J].Signal Processing, 
%                        2015,109: 165-171."

function w_LCSS=LCSS(RE,M,InterRange1,InterRange2,nullwidth,eta,a_s_n,SenLoc)

% INPUT ARGUEMENTS
% 'RE' is the sample covariance matrix of received signal
% 'M' is the number of elements
% 'InterRange1' is the first beam range of the interference
% 'InterRange2' is the second beam range of the interference
% 'nullwidth' is the width of the null
% 'eta' is the threshold of the output power in the null's sector
% 'a_s_n' is the nominal steering vector of desired signal

% OUTPUT ARGUEMENTS
% 'w_LCSS' is the estimate weight vector by LCSS method

dl=0.5;
derad=pi/180;
Q=zeros(M);
for theta_null=[InterRange1(1):0.1:InterRange1(2),InterRange2(1):0.1:InterRange2(2)]
    a_null=exp(-1j*2*pi*kron(SenLoc,sin(theta_null*derad)));
%     a_null=exp(-1j*(0:M-1)'*2*pi*dl*cos(theta_null*derad));
    Q=Q+a_null*a_null'*0.1*derad;
end
% Q=Q/(nullwidth*derad); % 公式（3）

[V,D]=eig(Q);
[d,ind]=sort(diag(D),'descend'); % sort the eigenvalues in descending order,d is eigenvalues in descending order,ind is a permutation vector of indices
Ds=D(ind,ind); % reorder the diagonal elements of D in descending order
Vs=V(:,ind); % reorder the columns of V to correspond to D

% 求使得q(r)<=eta成立的最小r值
for r=1:M
    if(r==1||q>eta)
         Vs_r=Vs(:,1:r); % 保留r个主特征值对应的特征向量
         P_ortho=eye(M)-RE^(-1/2)*Vs_r*inv(Vs_r'*inv(RE)*Vs_r)*Vs_r'*RE^(-1/2);
         q=0;
         for n=r+1:M
             q=q+Ds(n,n)*abs(a_s_n'*RE^(-1/2)*P_ortho*RE^(-1/2)*Vs(:,n))^2/...
                 (norm(P_ortho*RE^(-1/2)*a_s_n)^2);
         end
         r_opt=r; % 使得q(r)<=eta成立的最小r值
    end   
end

R_LCSSinv=RE^(-1/2)*P_ortho*RE^(-1/2); % 把Rin^(-1/2)*P_ortho*Rin^(-1/2)作为一个整体，...
                                      %        命名为R_LCSS
w_LCSS=R_LCSSinv*a_s_n/(a_s_n'*R_LCSSinv*a_s_n);
end