clc
clear
close all
SampFreq = 1000;
t = 0:1/SampFreq:1;

Sig1 = exp(j*(2*pi*(0*t)));


Sig2 = exp(j*2*pi*(150*t -300*t.^2));


Sig3 = exp(j*2*pi*(450*t -300*t.^2));

Sig4 = exp(j*2*pi*(-150*t +300*t.^2));

Sig5 = exp(j*2*pi*(-450*t + 300*t.^2));

Sig6 = exp(j*2*pi*(300*t -300*t.^2));

Sig7 = exp(j*2*pi*(-300*t + 300*t.^2));

Sig = Sig1 + Sig2 + Sig3 + Sig4 + Sig5 + Sig6 + Sig7;


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
bw = SampFreq/100;% use the Fourier model to smooth the detected ridge curves；orderIF1 could be larger than the following orderIF
orderIF = 10;
num = 7; % the number of the components
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
bw = SampFreq/50;
orderIF = 3;
orderamp = round(bw*length(Sig)/SampFreq);%瞬时幅值傅立叶级数阶次
[extr_Sig,ampmatrix,IFfit] = ICCD(Sig,SampFreq,f(findex),orderIF,orderamp,alpha);

%% finally fitted IFs
figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,IFfit,'linewidth',2.5);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
legend('C1','C2','C3','C4','C5','C6','C7')
set(gca,'FontSize',24)
set(gca,'linewidth',2);
axis([0 1 -SampFreq/2 SampFreq/2]);


figure
set(gcf,'Position',[20 100 640 200]);    
set(gcf,'Color','w'); 
plot(t,ampmatrix(1,:),'r','linewidth',4);
hold on
plot(t,real(extr_Sig(1,:)),'linewidth',2);
hold on
plot(t,real(Sig3 - extr_Sig(1,:)),'k','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('C1','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);

figure
set(gcf,'Position',[20 100 640 200]);    
set(gcf,'Color','w'); 
plot(t,ampmatrix(2,:),'r','linewidth',4);
hold on
plot(t,real(extr_Sig(2,:)),'linewidth',2);
hold on
plot(t,real(Sig4 - extr_Sig(2,:)),'k','linewidth',2);
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

figure
set(gcf,'Position',[20 100 640 200]);    
set(gcf,'Color','w'); 
plot(t,ampmatrix(4,:),'r','linewidth',4);
hold on
plot(t,real(extr_Sig(4,:)),'linewidth',2);
hold on
plot(t,real(Sig7 - extr_Sig(4,:)),'k','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('C4','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);

figure
set(gcf,'Position',[20 100 640 200]);    
set(gcf,'Color','w'); 
plot(t,ampmatrix(5,:),'r','linewidth',4);
hold on
plot(t,real(extr_Sig(5,:)),'linewidth',2);
hold on
plot(t,real(Sig5 - extr_Sig(5,:)),'k','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('C5','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);

figure
set(gcf,'Position',[20 100 640 200]);    
set(gcf,'Color','w'); 
plot(t,ampmatrix(6,:),'r','linewidth',4); %estimated Amplitude
hold on
plot(t,real(extr_Sig(6,:)),'linewidth',2);
hold on
plot(t,real(Sig6 - extr_Sig(6,:)),'k','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('C6','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);

figure
set(gcf,'Position',[20 100 640 200]);    
set(gcf,'Color','w'); 
plot(t,real(extr_Sig(7,:)),'linewidth',2);
hold on
plot(t,real(Sig1 - extr_Sig(7,:)),'k','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('C7','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);


