% Robust Beamforming Using PDL(projection transformation and diagnol loading)
% copyright Guo Yunzhou; Aug 25, 2020
% email: mailguoyunzhou@163.com

% according to paper "李文兴,毛晓军,孙亚秀.一种新的波束形成零陷展宽算法[J].
%                       电子与信息学报,2014,36(12):2882-2888."


function [R_PDL,w_PDL]=PDL_Li(M,InterRange1,InterRange2,K,RE,lambda,a_s_n,SenLoc)

% INPUT ARGUEMENTS
% 'M' is the number of elements
% 'InterRange1' is the first beam range of the interference
% 'InterRange2' is the second beam range of the interference
% 'K' is the number of base vectors
% 'RE' is the sample covariance matrix of received signal
% 'lambda' is the diagnol loading factor
% 'a_s_n' is the nominal steering vector of desired signal

% OUTPUT ARGUEMENTS
% 'R_PDL' is the estimate covariance matrix by PDL method
% 'w_PDL' is the estimate weight vector by PDL method

dl=0.5;
derad=pi/180;
R_theta=zeros(M);
for theta_proj=[InterRange1(1):0.1:InterRange1(2),InterRange2(1):0.1:InterRange2(2)]
    a_proj=exp(-1j*2*pi*kron(SenLoc,sin(theta_proj*derad)));
%     a_proj=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta_proj*derad));
    delta_theta=0.1; % 积分定义里的delta(x)
    R_theta=R_theta+a_proj*a_proj'*delta_theta;
end

[V,D]=eig(R_theta);
[d,ind]=sort(diag(D),'descend'); % sort the eigenvalues in descending order,d is eigenvalues in descending order,ind is a permutation vector of indices
Ds=D(ind,ind); % reorder the diagonal elements of D in descending order
Vs=V(:,ind); % reorder the columns of V to correspond to D
Ds=Ds(1:K,1:K); % maintain the signal-plus-interference subspace eigenvalues and remove the noise subspace eigenvalues
Vs=Vs(:,1:K); % maintain the signal-plus-interference subspace eigenvectors and remove the noise subspace eigenvectors
T=Vs*Vs'; % 投影算子

R_PDL=T*RE*T';    
R_PDL=R_PDL+lambda*eye(M);
w_PDL=inv(R_PDL)*a_s_n/(a_s_n'*inv(R_PDL)*a_s_n);
end