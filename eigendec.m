%Robust Adaptive Beamforming using WCPO method
%copyright Yunzhou Guo;Oct.11,2019
%email:690687612@qq.com
%
%[Vs,Ds]=eigendec(RE,ninter)
%
%INPUT ARGUMENTS
%RE is sample cov. matrix of received signal
%ninter is the number of interferences
%
%OUTPUT ARGUMENTS
%Vs is the signal-plus-interference subspace eigenvectors of RE
%Ds is the signal-plus-interference subspace eigenvalue matrix of RE


function w_EIG=eigendec(RE,ninter,a_s_n)
[V,D]=eig(RE);  %eigenvalue decomposition,V is the right eigenvector of RE,D is the eigenvalue matrix
[d,ind]=sort(diag(D),'descend'); %sort the eigenvalues in descending order,d is eigenvalues in descending order,ind is a permutation vector of indices
Ds=D(ind,ind); %reorder the diagonal elements of D in descending order
Ds=Ds(1:ninter+1,1:ninter+1); %maintain the signal-plus-interference subspace eigenvalues and remove the noise subspace eigenvalues
Vs=V(:,ind); %reorder the columns of V to correspond to D
Vs=Vs(:,1:ninter+1); %maintain the signal-plus-interference subspace eigenvectors and remove the noise subspace eigenvectors
w_EIG=Vs*inv(Ds)*Vs'*a_s_n;
end
        
        