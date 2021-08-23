% [R_CMT]=CMT(M,delta,RE)
%
% INPUT ARGUMENTS
% 'M' is the number of elements
% 'delta' is a parameter of CMT method
% 'RE' is the sample cov.
%
% OUTPUT ARGUMENTS
% 'R_CMT' is the modification of RE by CMT method
%
% Null broadening adaptive beamforming using CMT method
% copyright Guo Y; Feb 12,2020
% email:zhou78zhou66@163.com

function [R_CMT]=CMT(M,delta,RE)
for p=1:M
    for q=1:M
        T_CMT(p,q)=sinc((p-q)*delta);
    end
end
R_CMT=RE.*T_CMT;
