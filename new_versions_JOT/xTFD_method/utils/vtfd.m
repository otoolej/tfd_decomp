%-------------------------------------------------------------------------------
% vtfd: plot TFD and time/frequency plots
% 
% <<<<
% extracted from the memory-efficient fast TFD package:
% https://github.com/otoolej/memeff_TFDs
% (commit 87adc0dd09a25d941667c73c031713a02cbad533)
% >>>>
%
% Syntax: vtfd(tfd,s1,FS,n,k)
%
% Inputs:
%        tfd    - the time-frequency distribution 
%        s1     - time-domain signal
%        Fs     - sampling frequency
%        n      - time samples
%        k      - frequency samples
%
%

% John M. O' Toole, University College Cork
% Started: 27-08-2021
%
% last update: Time-stamp: <2021-08-27 16:48:01 (otoolej)>
%-------------------------------------------------------------------------------
function vtfd(tfd,s1,FS,n,k)
if(nargin<3 || isempty(FS)), FS=1; end
if(nargin<2 || isempty(s1)),
    TFD_ONLY=1;
    Ntime=size(tfd,1);  
else 
  TFD_ONLY=0; 
  Ntime=length(s1);
end
if(nargin<4 || isempty(n)), n=[]; end
if(nargin<5 || isempty(k)), k=[]; end


[N,M]=size(tfd);


% $$$ clf;


%---------------------------------------------------------------------
% Frame sizes:
%---------------------------------------------------------------------
X_AXIS_WIDTH=0.75;
Y_AXIS_HEIGHT=0.75;
X_AXIS_START=0.23;
Y_AXIS_START=0.23;
TIME_FREQ_PLOTS_WIDTH=0.15;
TIME_FREQ_PLOTS_GAP=0.07;


%---------------------------------------------------------------------
% X/Y tick labels
%---------------------------------------------------------------------
if(isempty(n))
    ntime=1:Ntime; ntime=ntime./FS;
    n=linspace(ntime(1),ntime(end),N);
else
    ntime=n;
end
Mh_time=ceil(Ntime/2);  Mh=ceil(M/2);
if(isempty(k))

    k=linspace(0,0.5,Mh);
    k_time=linspace(0,0.5,Mh_time);
    k=k.*FS;
else
    k_time=k;
end
% $$$ keyboard;

if(TFD_ONLY)
    imagesc(k,n,tfd); axis('xy');
    
else
  %---------------------------------------------------------------------
  % 1. time plot
  %---------------------------------------------------------------------
  s1=real(s1); 
  h_time=subplot(2,2,1);
  set(h_time,'position',[(X_AXIS_START-TIME_FREQ_PLOTS_WIDTH-TIME_FREQ_PLOTS_GAP) ...
                      Y_AXIS_START TIME_FREQ_PLOTS_WIDTH ...
                      Y_AXIS_HEIGHT]);
  plot(s1,1:length(s1)); 
  axis('tight');
  grid('on');
  set(h_time,'xticklabel',[]); set(h_time,'yticklabel',[]);


  %---------------------------------------------------------------------
  % 2. freq plot
  %---------------------------------------------------------------------
  S1=fft(s1);
  h_freq=subplot(2,2,4);
  set(h_freq,'position',[X_AXIS_START (Y_AXIS_START- ...
                                       TIME_FREQ_PLOTS_WIDTH-TIME_FREQ_PLOTS_GAP) ...
                      X_AXIS_WIDTH TIME_FREQ_PLOTS_WIDTH]);
  plot(1:Mh_time,abs(S1(1:Mh_time)).^2);
  axis('tight');
  grid('on');
  set(h_freq,'xticklabel',[]); set(h_freq,'yticklabel',[]);


  %---------------------------------------------------------------------
  % 2. time-freq plot
  %---------------------------------------------------------------------
  h_image=subplot(2,2,2);
  set(h_image,'position',[X_AXIS_START Y_AXIS_START X_AXIS_WIDTH ...
                      Y_AXIS_HEIGHT]);

  imagesc(k,n,tfd); axis('xy');
end

