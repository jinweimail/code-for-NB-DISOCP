clc;
clear all;
close all;

derad=pi/180;
M=10;
dl=0.5;

theta_s_r=0;
theta_s_n=0;
theta_i=[-40 50];

SenLoc = [0.1,0.4,0.9,1.6,2.1,2.4,2.9,3.6,4.1,4.6]';
a_s_r = exp(-1j*2*pi*kron(SenLoc,sin(theta_s_r*derad)));
a_s_n = exp(-1j*2*pi*kron(SenLoc,sin(theta_s_n*derad)));
a_i = exp(-1j*2*pi*kron(SenLoc,sin(theta_i*derad)));

Pn=1; %the noise power
PndB=10*log10(Pn);
PidB=[30 30]; %the interference to noise ratio
Pi=10.^(PidB/10);
Rin=a_i*(diag(Pi))*a_i'+Pn*eye(M);

%************* RCB参数 **************
epsilon_RCB=3;
% 
% %************* SQP参数 ***************
% BeamRange=[-5 5];
% theta1=BeamRange(1):0.0002:BeamRange(2);
% c1=exp(-1j*2*pi*kron(SenLoc,sin(theta1*derad)));
% % c1=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta1*derad));
% C=c1*c1'*((max(theta1*derad)-min(theta1*derad))/(length(theta1)-1)); % integral
% theta2=-90:0.0002:BeamRange(1);
% theta3=BeamRange(2):0.0002:90;
% c2=exp(-1j*2*pi*kron(SenLoc,sin(theta2*derad)));
% c3=exp(-1j*2*pi*kron(SenLoc,sin(theta3*derad)));
% % c2=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta2*derad));
% % c3=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta3*derad));
% C_bar=c2*c2'*((max(theta2*derad)-min(theta2*derad))/(length(theta2)-1))+...
%       c3*c3'*((max(theta3*derad)-min(theta3*derad))/(length(theta3)-1));
% sidelobe_level=abs(a_s_n'* C_bar* a_s_n);
% [V,D]=eig(C);
% [d,ind]=sort(diag(D),'descend');
% Vs=V(:,ind); % reorder the columns of V,the eigenvectors
% K_SQP=6; % the number of dominant eigenvalues of C
% U=Vs(:,1:K_SQP); % formula (11)
% Pa=eye(size(U*U'))-U*U'; % Pa is a projection matrix onto a subspace that is
%                          %    orthogonal to the actual steering vector
                         
%************* PDL参数 ***************
InterRange1=[-45 -35];
InterRange2=[45 55];
K_PDL=6; % 基向量个数
lambda_PDL=0.01;

%************* LCSS参数 **************
nullwidth=10; 
eta=1e-6;
InterRange1=[-45 -35];
InterRange2=[45 55];

%************* RAB-DISOCP-WC参数 ************
BeamRange=[-5 5];
Ripple=0.3;
StepSize = 0.5;
epsilon=0.1;
% threhold1 = 1e-5;

nsnapshot=500;
SNR=-15:5:15;
Nmc=2;
for iNmc=1:Nmc
    for iSNR=1:length(SNR) 
        PsdB=SNR(iSNR)+PndB;
        Ps=10^(PsdB/10);
        ss=sqrt(Ps)*sqrt(0.5)*(randn(1,nsnapshot)+j*randn(1,nsnapshot)); %SOI
        si=sqrt(diag(Pi))*sqrt(0.5)*(randn(2,nsnapshot)+j*randn(2,nsnapshot)); %the interference signal
        nE=sqrt(Pn)*sqrt(0.5)*(randn(M,nsnapshot)+j*randn(M,nsnapshot)); %noise
        xin=a_i*si+nE; %snapshot without signal
        xs=a_s_r*ss;
        x=xs+xin; %snapshot with signal
        REs=xs*xs'/nsnapshot;
        REin=xin*xin'/nsnapshot; %sample R of interference plus noise
        RE=x*x'/nsnapshot; %sample R of received signal
    
        SINR_OPT(iNmc,iSNR)=Ps*a_s_r'*inv(Rin)*a_s_r;
    
        w_SMI=inv(RE)*a_s_n/(a_s_n'*inv(RE)*a_s_n);
        SINR_SMI(iNmc,iSNR)=Ps*(abs(w_SMI'*a_s_r))^2/(w_SMI'*Rin*w_SMI);
%         SINR_SMI(iNmc,iSNR)=w_SMI'*REs*w_SMI/(w_SMI'*Rin*w_SMI);
        
        [a_RCB lambda_RCB] = RCBsphere(RE,a_s_n,epsilon_RCB);
        w_RCB=inv(RE)*a_RCB/(a_RCB'*inv(RE)*a_RCB);
        SINR_RCB(iNmc,iSNR)=Ps*(abs(w_RCB'*a_s_r))^2/(w_RCB'*Rin*w_RCB);
%         SINR_RCB(iNmc,iSNR)=w_RCB'*REs*w_RCB/(w_RCB'*Rin*w_RCB);
        
%         [w_SQP] =SQP(a_s_n, RE, Pa, C_bar, sidelobe_level);
%         SINR_SQP(iNmc,iSNR)=Ps*(abs(w_SQP'*a_s_r))^2/(w_SQP'*Rin*w_SQP);
%         SINR_SQP(iNmc,iSNR)=w_SQP'*REs*w_SQP/(w_SQP'*Rin*w_SQP);
       
        [R_PDL,w_PDL]=PDL_Mao(M,InterRange1,InterRange2,nullwidth,K_PDL,RE,lambda_PDL,a_s_n,SenLoc);
        SINR_PDL(iNmc,iSNR)=Ps*(abs(w_PDL'*a_s_r))^2/(w_PDL'*Rin*w_PDL);
%          SINR_PDL(iNmc,iSNR)=w_PDL'*REs*w_PDL/(w_PDL'*Rin*w_PDL);
        
        w_LCSS=LCSS(RE,M,InterRange1,InterRange2,nullwidth,eta,a_s_n,SenLoc);
        SINR_LCSS(iNmc,iSNR)=Ps*(abs(w_LCSS'*a_s_r))^2/(w_LCSS'*Rin*w_LCSS);
%         SINR_LCSS(iNmc,iSNR)=w_LCSS'*REs*w_LCSS/(w_LCSS'*Rin*w_LCSS);
        
        eta=1e-6;
        w_NB_DISOCP_WC = NB_DISOCP_WC_LA(RE,BeamRange,SenLoc,Ripple,StepSize,epsilon,InterRange1,InterRange2,eta);
        SINR_NB_DISOCP_WC(iNmc,iSNR)=Ps*(abs(w_NB_DISOCP_WC(:,end)'*a_s_r))^2/(w_NB_DISOCP_WC(:,end)'*Rin*w_NB_DISOCP_WC(:,end));
%         SINR_NB_DISOCP_WC(iNmc,iSNR)=w_NB_DISOCP_WC(:,end)'*REs*w_NB_DISOCP_WC(:,end)/(w_NB_DISOCP_WC(:,end)'*Rin*w_NB_DISOCP_WC(:,end));
    end
end

figure(1);
plot(SNR,10*log10(sum(SINR_OPT)/Nmc),'k-');hold on
plot(SNR,10*log10(sum(SINR_SMI)/Nmc),'k--');hold on
plot(SNR,10*log10(sum(SINR_RCB)/Nmc),'p--');hold on
% plot(SNR,10*log10(sum(SINR_SQP)/Nmc),'md--');hold on
plot(SNR,10*log10(sum(SINR_PDL)/Nmc),'gs-');hold on
plot(SNR,10*log10(sum(SINR_LCSS)/Nmc),'b^-');hold on
plot(SNR,10*log10(sum(SINR_NB_DISOCP_WC)/Nmc),'ro-');
legend('OPTIMAL SINR','SMI','RCB','PDL','LCSS','NB_DISOCP_WC');
xlabel('SNR（dB）');
ylabel('SINR(dB)');
axis([-15 15 -15 25]);
grid on;