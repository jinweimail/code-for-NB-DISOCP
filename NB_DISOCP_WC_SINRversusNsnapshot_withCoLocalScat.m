clc;
clear all;
close all;

derad=pi/180;
M=10;
dl=0.5;

theta_s_r=4;
% theta_s_n=0;
theta_i=[-40 50];

SenLoc = dl*[0:M-1]';
a_s_r = exp(-1j*2*pi*kron(SenLoc,sin(theta_s_r*pi/180)));
% a_s_n = exp(-1j*2*pi*kron(SenLoc,sin(theta_s_n*pi/180)));
a_i = exp(-1j*2*pi*kron(SenLoc,sin(theta_i*pi/180)));

Pn=1; %the noise power
PndB=10*log10(Pn);
PidB=[30 30]; %the interference to noise ratio
Pi=10.^(PidB/10);
SNR=15;
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

nsnapshot=[2:2:10,20:10:100];
Nmc=2;
for iNmc=1:Nmc
    a_s_n=a_s_r;
    for i_cls=1:4;
    theta_cls=normrnd(theta_s_r,8);
%     theta_cls=3-2*sqrt(3)+4*sqrt(3)*rand;
    a_cls=exp(-1j*2*pi*kron(SenLoc,sin(theta_cls*derad)));
%     a_cls=exp(-1j*(0:M-1)'*2*pi*dl*sin(theta_cls*derad));
    phi=2*pi*rand(1);
    a_s_n=a_s_n+exp(1j*phi)*a_cls;
    end
    for insnap=1:length(nsnapshot)
        ss=sqrt(Ps)*sqrt(0.5)*(randn(1,nsnapshot(insnap))+j*randn(1,nsnapshot(insnap))); %SOI
        si=sqrt(diag(Pi))*sqrt(0.5)*(randn(2,nsnapshot(insnap))+j*randn(2,nsnapshot(insnap))); %the interference signal
        nE=sqrt(Pn)*sqrt(0.5)*(randn(M,nsnapshot(insnap))+j*randn(M,nsnapshot(insnap))); %noise
        xin=a_i*si+nE; %snapshot without signal
        xs=a_s_r*ss;
        x=xs+xin; %snapshot with signal
        REs=xs*xs'/nsnapshot(insnap);
        REin=xin*xin'/nsnapshot(insnap); %sample R of interference plus noise
        RE=x*x'/nsnapshot(insnap); %sample R of received signal
    
        SINR_OPT(iNmc,insnap)=Ps*a_s_r'*inv(Rin)*a_s_r;
    
        w_SMI=inv(RE)*a_s_n/(a_s_n'*inv(RE)*a_s_n);
        SINR_SMI(iNmc,insnap)=Ps*(abs(w_SMI'*a_s_r))^2/(w_SMI'*Rin*w_SMI);
        
        w_EIG=eigendec(RE,ninter,a_s_n);
        SINR_EIG(iNmc,insnap)=Ps*(abs(w_EIG'*a_s_r))^2/(w_EIG'*Rin*w_EIG);
        
        [a_RCB lambda_RCB] = RCBsphere(RE,a_s_n,epsilon_RCB);
        w_RCB=inv(RE)*a_RCB/(a_RCB'*inv(RE)*a_RCB);
        SINR_RCB(iNmc,insnap)=Ps*(abs(w_RCB'*a_s_r))^2/(w_RCB'*Rin*w_RCB);

        [R_CMT]=CMT(M,delta_CMT,RE);
        w_CMT=inv(R_CMT)*a_s_n/(a_s_n'*inv(R_CMT)*a_s_n);
        SINR_CMT(iNmc,insnap)=Ps*(abs(w_CMT'*a_s_r))^2/(w_CMT'*Rin*w_CMT);
       
        [R_PDL,w_PDL]=PDL_Mao(M,InterRange1,InterRange2,nullwidth,K_PDL,RE,lambda_PDL,a_s_n,SenLoc);
        SINR_PDL(iNmc,insnap)=Ps*(abs(w_PDL'*a_s_r))^2/(w_PDL'*Rin*w_PDL);

        w_LCSS=LCSS(RE,M,InterRange1,InterRange2,nullwidth,eta_LCSS,a_s_n,SenLoc);
        SINR_LCSS(iNmc,insnap)=Ps*(abs(w_LCSS'*a_s_r))^2/(w_LCSS'*Rin*w_LCSS);
        
        w_NB_DISOCP_WC = NB_DISOCP_WC_LA(RE,BeamRange,SenLoc,Ripple,StepSize,epsilon,InterRange1,InterRange2,eta_NB_DISOCP);
        SINR_NB_DISOCP_WC(iNmc,insnap)=Ps*(abs(w_NB_DISOCP_WC(:,end)'*a_s_r))^2/(w_NB_DISOCP_WC(:,end)'*Rin*w_NB_DISOCP_WC(:,end));
    end
end

figure(1);
plot(nsnapshot,10*log10(sum(SINR_OPT)/Nmc),'k-','LineWidth',1);hold on
plot(nsnapshot,10*log10(sum(SINR_SMI)/Nmc),'+--','LineWidth',1);hold on
plot(nsnapshot,10*log10(sum(SINR_EIG)/Nmc),'x--','LineWidth',1);hold on
plot(nsnapshot,10*log10(sum(SINR_RCB)/Nmc),'p--','LineWidth',1);hold on
plot(nsnapshot,10*log10(sum(SINR_CMT)/Nmc),'md-','LineWidth',1);hold on
plot(nsnapshot,10*log10(sum(SINR_PDL)/Nmc),'gs-','LineWidth',1);hold on
plot(nsnapshot,10*log10(sum(SINR_LCSS)/Nmc),'b^-','LineWidth',1);hold on
plot(nsnapshot,10*log10(sum(SINR_NB_DISOCP_WC)/Nmc),'ro-','LineWidth',1);
legend('OPTIMAL SINR','SMI','EIG','RCB','CMT','PDL','LCSS','NB-DISOCP');
xlabel('Number of snapshot');
ylabel('SINR(dB)');
axis([0 100 -15 25]);
grid on;