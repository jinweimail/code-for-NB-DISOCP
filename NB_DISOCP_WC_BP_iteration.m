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

%************* NB-DISOCP参数 ************
BeamRange=[-5 5];
SenLoc = dl*[0:M-1]';
Ripple=0.3;
StepSize = 0.5;
epsilon=0.1;
eta_NB_DISOCP=1e-7;
InterRange1=[-45 -35];
InterRange2=[45 55];
% threhold1 = 1e-3;

nsnapshot=500;
theta=-90:90;
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
    
    w_NB_DISOCP_WC = NB_DISOCP_WC_LA(RE,BeamRange,SenLoc,Ripple,StepSize,epsilon,InterRange1,InterRange2,eta_NB_DISOCP);
    out_BP_NB_DISOCP_WC2(iNmc,:)=abs(w_NB_DISOCP_WC(:,2)'*A);
    out_BP_NB_DISOCP_WC3(iNmc,:)=abs(w_NB_DISOCP_WC(:,3)'*A);
    out_BP_NB_DISOCP_WC5(iNmc,:)=abs(w_NB_DISOCP_WC(:,5)'*A);
    out_BP_NB_DISOCP_WC11(iNmc,:)=abs(w_NB_DISOCP_WC(:,11)'*A);
    out_BP_NB_DISOCP_WC16(iNmc,:)=abs(w_NB_DISOCP_WC(:,16)'*A);
    out_BP_NB_DISOCP_WC20(iNmc,:)=abs(w_NB_DISOCP_WC(:,20)'*A);
    
    nIter=16;
    for inIter=1:nIter
        out_SINR_NB_DISOCP_WC(iNmc,inIter)=Ps*(abs(w_NB_DISOCP_WC(:,inIter)'*a_s_r))^2/(w_NB_DISOCP_WC(:,inIter)'*Rin*w_NB_DISOCP_WC(:,inIter));
    end
end

BP_NB_DISOCP_WC2=mean(out_BP_NB_DISOCP_WC2);
BP_NB_DISOCP_WC5=mean(out_BP_NB_DISOCP_WC5);
BP_NB_DISOCP_WC11=mean(out_BP_NB_DISOCP_WC11);
BP_NB_DISOCP_WC16=mean(out_BP_NB_DISOCP_WC16);

SINR_NB_DISOCP_WC=mean(out_SINR_NB_DISOCP_WC);

save('E:\研究\仿真\NB_DISOCP_WC\20200903 for paper')

figure(1);
plot(theta,20*log10(BP_NB_DISOCP_WC2/max(BP_NB_DISOCP_WC2)),'k:','LineWidth',1);hold on
plot(theta,20*log10(BP_NB_DISOCP_WC5/max(BP_NB_DISOCP_WC5)),'b--','LineWidth',1);hold on
plot(theta,20*log10(BP_NB_DISOCP_WC11/max(BP_NB_DISOCP_WC11)),'gs-','LineWidth',1);hold on
plot(theta,20*log10(BP_NB_DISOCP_WC16/max(BP_NB_DISOCP_WC16)),'r-','LineWidth',2);
legend('One Iteration','Four Iteration','Ten Iteration','Fifteen Iteration');
xlabel('Angle（degree）');
ylabel('Beampattern(dB)');
axis([-90 90 -90 0]);
% grid on;

figure(2);
plot(1:nIter-1,10*log10(SINR_NB_DISOCP_WC(2:nIter)),'ro-','LineWidth',2);
xlabel('Iteration number');
ylabel('OUTPUT SINR(dB)');
% grid on