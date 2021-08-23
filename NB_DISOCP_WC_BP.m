clc;
clear all;
close all;

derad=pi/180;
M=10;
dl=0.5;

theta_s_r=0;
theta_s_n=0;
theta_i=[-40 50];

a_s_r=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta_s_r*derad));
a_s_n=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta_s_n*derad));
a_i=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta_i*derad));

Pn=1; %the noise power
PndB=10*log10(Pn);
PidB=[30 30]; %the interference to noise ratio
Pi=10.^(PidB/10);
SNR=0;
PsdB=SNR+PndB;
Ps=10^(PsdB/10);
Rin=a_i*(diag(Pi))*a_i'+Pn*eye(M);

%************* EIG参数 **************
ninter=length(theta_i);

%************* RCB参数 **************
epsilon_RCB=3;

%************* CMT参数 ***************
delta_CMT=0.05;

%************* PDL参数 ***************
InterRange1=[-45 -35];
InterRange2=[45 55];
K_PDL=6; % 基向量个数
lambda_PDL=0.01;

%************* LCSS参数 **************
nullwidth=10; 
eta_LCSS=1e-7;
InterRange1=[-45 -35];
InterRange2=[45 55];

%************* NB-DISOCP参数 ************
BeamRange=[-5 5];
SenLoc = dl*[0:M-1]';
Ripple=0.3;
StepSize = 0.5;
epsilon=0.1;
eta_NB_DISOCP=1e-7;
% threhold1 = 1e-3;

nsnapshot=500;
theta=-90:0.1:90;
A=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta*derad));
Nmc=2;
for iNmc=1:Nmc
    ss=sqrt(Ps)*sqrt(0.5)*(randn(1,nsnapshot)+j*randn(1,nsnapshot)); %SOI
    si=sqrt(diag(Pi))*sqrt(0.5)*(randn(2,nsnapshot)+j*randn(2,nsnapshot)); %the interference signal
    nE=sqrt(Pn)*sqrt(0.5)*(randn(M,nsnapshot)+j*randn(M,nsnapshot)); %noise
    xin=a_i*si+nE; %snapshot without signal
    x=a_s_r*ss+xin; %snapshot with signal
    REin=xin*xin'/nsnapshot; %sample R of interference plus noise
    RE=x*x'/nsnapshot; %sample R of received signal
    
    w_SMI=inv(RE)*a_s_n/(a_s_n'*inv(RE)*a_s_n);
    BP_SMI(iNmc,:)=abs(w_SMI'*A);
    
    w_EIG=eigendec(RE,ninter,a_s_n);
    BP_EIG(iNmc,:)=abs(w_EIG'*A);
    
    [a_RCB lambda_RCB] = RCBsphere(RE,a_s_n,epsilon_RCB);
    w_RCB=inv(RE)*a_RCB/(a_RCB'*inv(RE)*a_RCB);
    BP_RCB(iNmc,:)=abs(w_RCB'*A);
    
    [R_CMT]=CMT(M,delta_CMT,RE);
    w_CMT=inv(R_CMT)*a_s_n/(a_s_n'*inv(R_CMT)*a_s_n);
    BP_CMT(iNmc,:)=abs(w_CMT'*A);
    
    [R_PDL,w_PDL]=PDL_Mao(M,InterRange1,InterRange2,nullwidth,K_PDL,RE,lambda_PDL,a_s_n,SenLoc);
    BP_PDL(iNmc,:)=abs(w_PDL'*A);
    
    w_LCSS=LCSS(RE,M,InterRange1,InterRange2,nullwidth,eta_LCSS,a_s_n,SenLoc);
    BP_LCSS(iNmc,:)=abs(w_LCSS'*A);
    
    w_NB_DISOCP_WC = NB_DISOCP_WC_LA(RE,BeamRange,SenLoc,Ripple,StepSize,epsilon,InterRange1,InterRange2,eta_NB_DISOCP);
    BP_NB_DISOCP_WC(iNmc,:)=abs(w_NB_DISOCP_WC(:,end)'*A);
end

figure(1);
% plot(theta,20*log10((sum(BP_SMI)/Nmc)/max(sum(BP_SMI)/Nmc)),'k:','LineWidth',1);hold on
% plot(theta,20*log10((sum(BP_EIG)/Nmc)/max(sum(BP_EIG)/Nmc)),'--');hold on
% plot(theta,20*log10((sum(BP_RCB)/Nmc)/max(sum(BP_RCB)/Nmc)),'.-');hold on
plot(theta,20*log10((sum(BP_CMT)/Nmc)/max(sum(BP_CMT)/Nmc)),'b:','LineWidth',2);hold on
plot(theta,20*log10((sum(BP_PDL)/Nmc)/max(sum(BP_PDL)/Nmc)),'g--','LineWidth',2);hold on
plot(theta,20*log10((sum(BP_LCSS)/Nmc)/max(sum(BP_LCSS)/Nmc)),'m-.','LineWidth',2);hold on
plot(theta,20*log10((sum(BP_NB_DISOCP_WC)/Nmc)/max(sum(BP_NB_DISOCP_WC)/Nmc)),'r-','LineWidth',3); 
% legend('SMI','EIG','RCB','CMT','PDL','LCSS','NB-DISOCP');
legend('CMT','PDL','LCSS','NB-DISOCP');
% xlabel('Angle（degree）');
% ylabel('Beampattern(dB)');
xlabel('角度（度）');
ylabel('增益(dB)');
axis([-90 90 -100 0]);
% grid on;