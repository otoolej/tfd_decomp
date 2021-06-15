clc
clear
close all
SampFreq = 1000;
t = 0:1/SampFreq:1;

Sig1 = exp(j*(2*pi*(0*t))); %


Sig2 = exp(-0.5*t).*exp(j*(140*pi*t + 150*sin(2*pi*t)));


Sig3 = 0.8*exp(0.5*t).*exp(j*(-140*pi*t - 150*sin(2*pi*t)));

Sig = Sig1 + Sig2 + Sig3;







figure
set(gcf,'Position',[20 100 320 250]);	    
set(gcf,'Color','w'); 
plot(t,Sig,'linewidth',1);
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',12)
set(gca,'linewidth',1);
axis([0 1 -1.2 3.2])

%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
Ratio = 0;
window = 128;
Nfrebin = 1024;
figure
[Spec,f] = STFT(Sig(:),SampFreq,Nfrebin,window);
imagesc(t,f,abs(Spec)); 
axis([0 1 -SampFreq/2 SampFreq/2]);
set(gcf,'Position',[20 100 320 250]);	 
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',12);
set(gcf,'Color','w');	
%% ridge extraction and fitting
bw = SampFreq/100;%% use the Fourier model to smooth the detected ridge curves£»orderIF1 could be larger than the following orderIF
orderIF = 20;
num = 3; % the number of the components
delta = 20;
alpha = 5;
[fidexmult, tfdv] = extridge_mult(Sig, SampFreq, num, delta, orderIF,bw,Nfrebin,window,alpha);


figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,f(fidexmult),'linewidth',3);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
axis([0 1 -SampFreq/2 SampFreq/2]);



%% ridge path regrouping
thrf = length(f)/30;

tic
[findex,interset] = RPRG(fidexmult,thrf);
toc


figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,f(findex),'linewidth',3);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
axis([0 1 -490 490]);


%% signal decomposition
alpha = 0.5;
bw = SampFreq/40; % bandwidth of the ICCD
orderamp = round(bw*length(Sig)/SampFreq);%% Fourier order for characterizing signal amplitudes
[extr_Sig,ampmatrix,IFfit] = ICCD(Sig,SampFreq,f(findex),orderIF,orderamp,alpha);% ICCD performs the joint-optimization scheme using all the obtained IFs

%% finally fitted IFs
figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,IFfit,'linewidth',3);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
axis([0 1 -SampFreq/2 SampFreq/2]);

%% reconstructed components
figure
set(gcf,'Position',[20 100 640 200]);    
set(gcf,'Color','w'); 
plot(t,real(extr_Sig(1,:)),'linewidth',2);
hold on
plot(t,real(Sig1 - extr_Sig(1,:)),'k','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('C1','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);

figure
set(gcf,'Position',[20 100 640 200]);
set(gcf,'Color','w'); 
plot(t,ampmatrix(2,:),'r','linewidth',4);% estimated amplitude
hold on
plot(t,real(extr_Sig(2,:)),'linewidth',2); % estimated signal
hold on
plot(t,real(Sig3 - extr_Sig(2,:)),'k','linewidth',2); % estimation errors
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('C2','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);

figure
set(gcf,'Position',[20 100 640 200]);   
set(gcf,'Color','w'); 
plot(t,ampmatrix(3,:),'r','linewidth',4);
hold on
plot(t,real(extr_Sig(3,:)),'linewidth',2);
hold on
plot(t,real(Sig2 - extr_Sig(3,:)),'k','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('C3','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);




