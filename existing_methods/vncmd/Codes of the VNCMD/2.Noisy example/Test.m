%%%%%%%%%%%%%%%%% signal with very strong background noise %%%%%%%%%%%%%%%%
clc
clear
close all
SampFreq = 500;
t = 0:1/SampFreq:6;
Sig1 = cos(2*pi*(1.3+25*t + 4*t.^2 - 0.8*t.^3 + 0.07*t.^4));
IF1 = 25 + 8*t - 2.4*t.^2 + 0.28*t.^3;
Sig2 = cos(2*pi*(2.6+40*t + 8*t.^2 -1.6*t.^3 + 0.14*t.^4));
IF2 = 40 + 16*t - 4.8*t.^2 + 0.56*t.^3;

Sig = Sig1 + Sig2;

STD = 1;% standard deviation of the added noise; 
noise = addnoise(length(Sig),0,STD); % the mean value of the noise is zero
Sig = Sig + noise; % the SNRs of the two noisy modes are around -3 dB: SNR(Sig1, Sig1+noise)

figure
set(gcf,'Position',[20 100 320 250]);	    
set(gcf,'Color','w'); 
plot(t,Sig,'linewidth',1);
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',12)
set(gca,'linewidth',1);
axis([0 6 -5 5])
%% STFT
Ratio = 0;
figure
[Spec,f] = STFT(Sig',SampFreq,512,256);
imagesc(t,f,abs(Spec)); 
axis([0 6 0 100]);
set(gcf,'Position',[20 100 320 250]);	 
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',12);
set(gcf,'Color','w');	
%% parameter setting
iniIF = [30*ones(1,length(t));60*ones(1,length(t))];% initial guess for the IFs for the three signal modes
alpha = 3e-4; %or alpha = 2e-4; % if this parameter is larger, it will help the algorithm to find correct modes even the initial IFs are too rough. But it will introduce more noise and also may increase the interference between the signal modes
beta = 1e-9; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
var = STD^2;% noise variance
tol = 1e-8;%

[IFmset IA smset] = VNCMD(Sig,SampFreq,iniIF,alpha,beta,var,tol);
%% Relative errors of the finally estimated IFs
   RE1 =  norm(IFmset(1,:,end)-IF1)/norm(IF1)
   RE2 =  norm(IFmset(2,:,end)-IF2)/norm(IF2)
%% estimated IF
figure
plot(t,[IF1;IF2],'b','linewidth',3) % true IFs
hold on
plot(t,IFmset(:,:,end),'r','linewidth',3) % finally estimated IFs
set(gcf,'Position',[20 100 640 500]);	 
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');	
axis([0 6 0 100])
%% Reconstructed modes
figure
set(gcf,'Position',[20 100 640 200]);	    
set(gcf,'Color','w'); 
plot(t,smset(1,:,end),'linewidth',2)	% estimated mode
hold on
plot(t,Sig1 - smset(1,:,end),'k','linewidth',2) % estimation errors
hold on
plot(t,IA(1,:),'r','linewidth',3) % estimated IAs
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('m1','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);	
axis([0 6 -1.5 1.5])

figure
set(gcf,'Position',[20 100 640 200]);	    
set(gcf,'Color','w'); 
plot(t,smset(2,:,end),'linewidth',2)	% estimated mode
hold on
plot(t,Sig2 - smset(2,:,end),'k','linewidth',2)  % estimation errors
hold on
plot(t,IA(2,:),'r','linewidth',3) % estimated IAs
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('m2','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);	
axis([0 6 -1.5 1.5])
%% SNRs of the reconstructed modes
snr1 = SNR(Sig1,smset(1,:,end))
snr2 = SNR(Sig2,smset(2,:,end))

Residue = Sig - sum(smset(:,:,end)); % Residual signals

stdres = std(Residue) % standard deviation of the residue; if the algorithm converges to correct results, this value will exactly equal to the standard deviation of the noise


