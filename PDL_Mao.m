function [R_PDL,w_PDL]=PDL_Mao(M,InterRange1,InterRange2,nullwidth,K,RE,lambda,a_s_n,SenLoc)

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

derad=pi/180;
theta_proj_mao=[InterRange1(1):1:InterRange1(2),InterRange2(1):1:InterRange2(2)];
    for itheta_proj_mao=1:length(theta_proj_mao)
        a_proj_mao=exp(-1j*2*pi*kron(SenLoc,sin(theta_proj_mao(itheta_proj_mao)*derad)));
%         a_proj_mao=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta_proj_mao(itheta_proj_mao)*derad));
        P_mao(itheta_proj_mao)=1./(a_proj_mao'*inv(RE)*a_proj_mao); % 计算干扰区间上的功率谱
    end
    [pks,locs] = findpeaks(abs(P_mao)); % 寻找干扰区间上功率谱函数的峰值，pks为峰值，locs为峰值所在位置
    theta_p_mao=theta_proj_mao(locs); % 寻找到的干扰位置
    
    z=zeros(M);
    for itheta_p_mao=1:length(theta_p_mao)
        a_theta_p_mao=exp(-1j*2*pi*kron(SenLoc,sin(theta_p_mao(itheta_p_mao)*derad)));
%         a_theta_p_mao=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta_p_mao(itheta_p_mao)*derad));
        z=z+a_theta_p_mao*a_theta_p_mao';
    end
    [R_theta_mao]=CMT(M,0.3*nullwidth*derad,z); % 见文献式（11）注
    [V_mao,D_mao]=eig(R_theta_mao);
    [d_mao,ind_mao]=sort(diag(D_mao),'descend'); % sort the eigenvalues in descending order,d is eigenvalues in descending order,ind is a permutation vector of indices
    Ds_mao=D_mao(ind_mao,ind_mao); % reorder the diagonal elements of D in descending order
    Vs_mao=V_mao(:,ind_mao); % reorder the columns of V to correspond to D
   
    Ds_mao=Ds_mao(1:K,1:K); % maintain the signal-plus-interference subspace eigenvalues and remove the noise subspace eigenvalues
    Vs_mao=Vs_mao(:,1:K); % maintain the signal-plus-interference subspace eigenvectors and remove the noise subspace eigenvectors
    T_mao=Vs_mao*Vs_mao'; % 投影算子
    R_proj_mao=T_mao*RE*T_mao';    
   
    R_PDL=R_proj_mao+lambda*eye(M);
    w_PDL=inv(R_PDL)*a_s_n/(a_s_n'*inv(R_PDL)*a_s_n);
end