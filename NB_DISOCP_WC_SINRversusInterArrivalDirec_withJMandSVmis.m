clc;
clear all;
close all;

derad=pi/180;
M=10;
dl=0.5;

theta_s_r=0;
theta_s_n=7;
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

nsnapshot=100;
inter_move=-7:7;
Nmc=2;
for iNmc=1:Nmc
    for inim=1:length(inter_move)
        ss=sqrt(Ps)*sqrt(0.5)*(randn(1,nsnapshot)+j*randn(1,nsnapshot)); %SOI
        si=sqrt(diag(Pi))*sqrt(0.5)*(randn(2,nsnapshot)+j*randn(2,nsnapshot)); %the interference signal
        nE=sqrt(Pn)*sqrt(0.5)*(randn(M,nsnapshot)+j*randn(M,nsnapshot)); %noise
        xin=a_i*si+nE; %snapshot without signal
        x=a_s_r*ss+xin; %snapshot with signal
        REin=xin*xin'/nsnapshot; %sample R of interference plus noise
        RE=x*x'/nsnapshot; %sample R of received signal
    
        theta_i_move=[-40+inter_move(inim) 50-inter_move(inim)]; % the DOA of interference
        a_i_move=exp(-j*(0:M-1)'*2*pi*dl*sin(theta_i_move*derad));
        Rin_move=a_i_move*diag(Pi)*a_i_move'+Pn*eye(M);
        xin_move=a_i_move*si+nE; %snapshot without signal
        x_move=a_s_r*ss+xin_move; %snapshot with signal
        REin_move=xin_move*xin_move'/nsnapshot; %sample R of interference plus noise
        RE_move=x_move*x_move'/nsnapshot;
    
        SINR_OPT(iNmc,inim)=Ps*a_s_r'*inv(Rin)*a_s_r;
    
        w_SMI=inv(RE)*a_s_n/(a_s_n'*inv(RE)*a_s_n);
        SINR_SMI(iNmc,inim)=Ps*(abs(w_SMI'*a_s_r))^2/(w_SMI'*Rin_move*w_SMI);
        
        w_EIG=eigendec(RE,ninter,a_s_n);
        SINR_EIG(iNmc,inim)=Ps*(abs(w_EIG'*a_s_r))^2/(w_EIG'*Rin_move*w_EIG);
        
        [a_RCB lambda_RCB] = RCBsphere(RE,a_s_n,epsilon_RCB);
        w_RCB=inv(RE)*a_RCB/(a_RCB'*inv(RE)*a_RCB);
        SINR_RCB(iNmc,inim)=Ps*(abs(w_RCB'*a_s_r))^2/(w_RCB'*Rin_move*w_RCB);
        
        [R_CMT]=CMT(M,delta_CMT,RE);
        w_CMT=inv(R_CMT)*a_s_n/(a_s_n'*inv(R_CMT)*a_s_n);
        SINR_CMT(iNmc,inim)=Ps*(abs(w_CMT'*a_s_r))^2/(w_CMT'*Rin_move*w_CMT);
       
        [R_PDL,w_PDL]=PDL_Mao(M,InterRange1,InterRange2,nullwidth,K_PDL,RE,lambda_PDL,a_s_n,SenLoc);
        SINR_PDL(iNmc,inim)=Ps*(abs(w_PDL'*a_s_r))^2/(w_PDL'*Rin_move*w_PDL);
        
        w_LCSS=LCSS(RE,M,InterRange1,InterRange2,nullwidth,eta_LCSS,a_s_n,SenLoc);
        SINR_LCSS(iNmc,inim)=Ps*(abs(w_LCSS'*a_s_r))^2/(w_LCSS'*Rin_move*w_LCSS);
       
        w_NB_DISOCP_WC = NB_DISOCP_WC_LA(RE,BeamRange,SenLoc,Ripple,StepSize,epsilon,InterRange1,InterRange2,eta_NB_DISOCP);
        SINR_NB_DISOCP_WC(iNmc,inim)=Ps*(abs(w_NB_DISOCP_WC(:,end)'*a_s_r))^2/(w_NB_DISOCP_WC(:,end)'*Rin_move*w_NB_DISOCP_WC(:,end));
end
end

figure(1);
plot(inter_move,10*log10(sum(SINR_OPT)/Nmc),'k-','LineWidth',1);hold on
plot(inter_move,10*log10(sum(SINR_SMI)/Nmc),'+--','LineWidth',1);hold on
plot(inter_move,10*log10(sum(SINR_EIG)/Nmc),'x--','LineWidth',1);hold on
plot(inter_move,10*log10(sum(SINR_RCB)/Nmc),'p--','LineWidth',1);hold on
plot(inter_move,10*log10(sum(SINR_CMT)/Nmc),'md-','LineWidth',1);hold on
plot(inter_move,10*log10(sum(SINR_PDL)/Nmc),'gs-','LineWidth',1);hold on
plot(inter_move,10*log10(sum(SINR_LCSS)/Nmc),'b^-','LineWidth',1);hold on
plot(inter_move,10*log10(sum(SINR_NB_DISOCP_WC)/Nmc),'ro-','LineWidth',1);    
legend('OPTIMAL SINR','SMI','EIG','RCB','CMT','PDL','LCSS','NB-DISOCP');
% xlabel('干扰角度偏移（度）');
% ylabel('增益(dB)');
xlabel('Deviation of interference arrival direction（degree）');
ylabel('OUTPUT SINR(dB)');
axis([-7 7 -45 12]);
grid on;