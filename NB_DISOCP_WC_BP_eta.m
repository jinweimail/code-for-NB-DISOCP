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

%************* NB-DISOCP²ÎÊý ************
BeamRange=[-5 5];
SenLoc = dl*[0:M-1]';
Ripple=0.3;
StepSize = 0.5;
epsilon=0.1;
eta=[1e-4,1e-5,1e-6,1e-7,1e-8,1e-9];
% eta=1e-4;
InterRange1=[-45 -35];
InterRange2=[45 55];
% threhold1 = 1e-4;

nsnapshot=100;
theta=-90:90;
A=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta*derad));


    ss=sqrt(Ps)*sqrt(0.5)*(randn(1,nsnapshot)+j*randn(1,nsnapshot)); %SOI
    si=sqrt(diag(Pi))*sqrt(0.5)*(randn(2,nsnapshot)+j*randn(2,nsnapshot)); %the interference signal
    nE=sqrt(Pn)*sqrt(0.5)*(randn(M,nsnapshot)+j*randn(M,nsnapshot)); %noise
    xin=a_i*si+nE; %snapshot without signal
    x=a_s_r*ss+xin; %snapshot with signal
    REin=xin*xin'/nsnapshot; %sample R of interference plus noise
    RE=x*x'/nsnapshot; %sample R of received signal
    
for ieta=1:length(eta)
    w_NB_DISOCP_WC = NB_DISOCP_WC_LA(RE,BeamRange,SenLoc,Ripple,StepSize,epsilon,InterRange1,InterRange2,eta(ieta));
    BP_NB_DISOCP_WC(ieta,:)=abs(w_NB_DISOCP_WC(:,end)'*A);
end


figure(1);
plot(theta,20*log10(BP_NB_DISOCP_WC(1,:)/max(BP_NB_DISOCP_WC(1,:))),'k:','LineWidth',1.1);hold on
plot(theta,20*log10(BP_NB_DISOCP_WC(2,:)/max(BP_NB_DISOCP_WC(2,:))),'b--','LineWidth',1.1);hold on
plot(theta,20*log10(BP_NB_DISOCP_WC(3,:)/max(BP_NB_DISOCP_WC(3,:))),'-.','LineWidth',1.1);hold on
plot(theta,20*log10(BP_NB_DISOCP_WC(4,:)/max(BP_NB_DISOCP_WC(4,:))),'m-','LineWidth',1.1);hold on
plot(theta,20*log10(BP_NB_DISOCP_WC(5,:)/max(BP_NB_DISOCP_WC(5,:))),'g.-','LineWidth',1.1);hold on
plot(theta,20*log10(BP_NB_DISOCP_WC(6,:)/max(BP_NB_DISOCP_WC(6,:))),'r-','LineWidth',1.1);
legend('\xi^2=10^{-4}','\xi^2=10^{-5}','\xi^2=10^{-6}','\xi^2=10^{-7}','\xi^2=10^{-8}','\xi^2=10^{-9}');
xlabel('Angle£¨degree£©');
ylabel('Beampattern(dB)');
axis([-90 90 -100 0]);
grid on;
