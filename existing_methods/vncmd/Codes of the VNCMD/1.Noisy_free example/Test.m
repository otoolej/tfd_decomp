%%%%%%%%%%%%%% signal with intersecting IFs %%%%%%%%%%%%%%%%%
clc
clear
close all
SampFreq = 2000;
t = 0:1/SampFreq:1;
Sig1 = (1+0.5*cos(2*pi*t)).*cos(2*pi*(0.2 + 532*t -474*t.^2 + 369*t.^3));
IF1 = 532 - 948*t + 1107*t.^2;
Sig2 = (1+0.5*cos(2*pi*t)).*cos(2*pi*(0.8+50*t + 525*t.^2 -300*t.^3));
IF2 = 50 + 1050*t - 900*t.^2;
Sig = Sig1+Sig2;
figure
set(gcf,'Position',[20 100 320 250]);	    
set(gcf,'Color','w'); 
plot(t,Sig,'linewidth',1);
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'FontSize',12)
set(gca,'linewidth',1);
axis([0 1 -4 4])
%% STFT
Ratio = 0;
figure
[Spec,f] = STFT(Sig',SampFreq,512,256);
imagesc(t,f,abs(Spec)); 
axis([0 1 0 700]);
set(gcf,'Position',[20 100 320 250]);	 
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'YDir','normal')
set(gca,'FontSize',12);
set(gcf,'Color','w');	
%% parameter setting
iniIF = [700*ones(1,length(t));20*ones(1,length(t))];% initial guess for the IFs for the three signal modes
alpha = 5e-6; 
beta = 1e-6; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
var = 0;% noise variance
tol = 1e-8;%
tic
[IFmset IA smset] = VNCMD(Sig,SampFreq,iniIF,alpha,beta,var,tol);
toc

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
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');	
axis([0 1 0 800])
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
axis([0 1 -1.9 1.9])

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
axis([0 1 -1.9 1.9])


