% Robust Beamforming Using RAB-DISOCP-WC  for linear array
% copyright Jin Wei; Nov 13, 2011
% email: jinzhengwei_521@163.com

%according to paper"Fan Wang, Venkataramanan Balakrishnan, Philip Yuanping Zhou, Jason Jiangnan Chen, Roy Yang, 
%                    and Colin Frank, Optimal Array Pattern Synthesis Using Semidefinite Programming,"IEEE Transactions on 
%                    Signal Processing, vol. 51, No.5, May 2003"



function w_NB_DISOCP_WC = NB_DISOCP_WC_LA(R,BeamRange,SenLoc,Ripple,StepSize,epsilon,InterRange1,InterRange2,eta)

% INPUT ARGUMENTS
% 'R'   the covariance matrix of received snapshots
% 'BeamRange'    is beam range of the mainlobe to control
% 'SenLoc'    is the position of the elements, eg:[0 0.5 1 1.5 2]
% 'Ripple'   is the ripple of the mainlobe to control
% 'Sampling'  is the samp number in [0,pi]
% 'epsilon' is the relative reguler factor
% 'threhold1' is the stopping threhold for the iteration process

% OUTPUT ARGUMENT
% 'w_RB_DISOCP_WC'  is the estimate weight vector by RB-DISOCP-WC method
%


[M,M] = size(R);

epsilon = real(epsilon*R(1,1));

theta = BeamRange(1):StepSize:BeamRange(2);
theta_null1 = InterRange1(1):StepSize:InterRange1(2);
theta_null2 = InterRange2(1):StepSize:InterRange2(2);


tao = kron(SenLoc,sin(theta*pi/180)) ;
A = exp(-1j*2*pi*tao);
A_null1=exp(-1j*2*pi*kron(SenLoc,sin(theta_null1*pi/180)));
A_null2=exp(-1j*2*pi*kron(SenLoc,sin(theta_null2*pi/180)));


% dBsidelobe = 0;
% dBmainlobe = -30;
U = 10^(+Ripple/20);
L = 10^(-Ripple/20);
N = size(A,2);
N_null1=size(A_null1,2);
N_null2=size(A_null2,2);

dl=0.5;
derad=pi/180;

%********************************************************************
% formulate and solve the magnitude design problem
%********************************************************************
cvx_begin
  cvx_quiet(1);
  variable w(M,1) complex
  if epsilon==0
      minimize ( quad_form(w,R))
  else
      minimize ( quad_form(w,R+epsilon*eye(M,M)))
  end
%   expression e(N)
%         for i=1:N
%            e(i)=norm(w'*A(:,i)-1,2)/N;
%         end
  subject to
%         sum(e) <= delta_max^2;
%         for i = 1:size(Asidelobe,2)
%             norm(w'*Asidelobe(:,i)) <= 10^(dBsidelobe/20);
%         end
%         for i = 1:N
%             norm(w'*A(:,i)-1) <= 10^(dBmainlobe/20);
% %             norm(w'*A(:,i)-1) <= 0.01;
%         end
      
        for i = 1:N
            norm(w'*A(:,i)-(U+L)/2) <= (U-L)/2;
        end
        for i_null1 = 1:N_null1
        quad_form(w,A_null1(:,i_null1)*A_null1(:,i_null1)')<=eta;
        end
        for i_null2 = 1:N_null2
        quad_form(w,A_null2(:,i_null2)*A_null2(:,i_null2)')<=eta;
        end
%         norm(w) <= 1;
cvx_end

w_NB_DISOCP_WC(:,1) = w;
count = 2;
% w1 = w/2;
while count <= 20
% while (count==2||norm(w1-w2) > threhold1)
    w1 = w/2;
    cvx_begin
      cvx_quiet(1);
      variable w2(M,1) complex
      if epsilon==0
          minimize ( quad_form(w1+w2,R))
      else
          minimize ( quad_form(w1+w2,R+epsilon*eye(M,M)))
      end
      subject to
        for i=1:N
            quad_form(w1+w2,A(:,i)*A(:,i)') <= U^2;
            4*real(w1'*A(:,i)*A(:,i)'*w2) >= L^2;
        end
        for i_null1 = 1:N_null1
        quad_form(w1+w2,A_null1(:,i_null1)*A_null1(:,i_null1)')<=eta;
        end
        for i_null2 = 1:N_null2
        quad_form(w1+w2,A_null2(:,i_null2)*A_null2(:,i_null2)')<=eta;
        end
    cvx_end
    w = w1+w2;
    w_NB_DISOCP_WC(:,count) = w;

%     if norm(w1-w2) <= threhold1
%         break;
%     end
%     w1 = w2;
%     w_MC(:,count) = w1+w2;
    count = count+1;
end



