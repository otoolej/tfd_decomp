%-------------------------------------------------------------------------------
% USE: sig=gen_signal(TYPE,M,A,L,K,E)
% 
%
%  TYPE is the type of signal, one of { 'PWLFM', 'AMP_VAR_PWLFM', 'APPROX_PWLFM' }
%  M is the number of pieces for PWLFM
%  A is vector of slopes for the number of pieces 
%    i.e. [A_1,A_2,...,A_M]
%  L is vector of lengths for the number of pieces 
%    i.e. [L_1,L_2,...,L_M]
%  K is the number of harmonic components
%  E is vector of amplitudes for the number of harmonic components
%    i.e. [E_1,E_2,...,E_M]
%  p is vector of phase-offsets for the each of harmonic component (and fundamental)
%    i.e. [p_1,p_2,...,p_{M+1}]
%  N number of samples 
%  epoch_length: epoch length (in seconds)
%  E_turning_points: if time-varying amplitude, these represent the location (in
%  seconds) for the turning points to constuct a time-varying amplidute function
%  E_tv_points: values for turning points 
%-------------------------------------------------------------------------------

%
%   Copyright (c) 2010, John M. O' Toole, The University of Queensland
%   All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following
%  conditions are met:
%      * Redistributions of source code must retain the above
%        copyright notice, this list of conditions and the following
%        disclaimer.
%      * Redistributions in binary form must reproduce the above
%        copyright notice, this list of conditions and the following
%        disclaimer in the documentation and/or other materials
%        provided with the distribution.
%      * Neither the name of the The University of Queensland nor the 
%        names of its contributors may be used to endorse or promote 
%        products derived from this software without specific prior 
%        written permission.
%  
%  THIS SOFTWARE IS PROVIDED BY JOHN M. O' TOOLE ''AS IS'' AND ANY
%  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JOHN M. O' TOOLE BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
%  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
%  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%  DAMAGE.
%
%-------------------------------------------------------------------------------
function [y,iflaw]=gen_signal(TYPE,M,A,L,K,E,p,N,epoch_length,TV_AMP,E_turning_points,E_tv_points,f0)

% set the defaults:
if(nargin<1 || isempty(TYPE)) TYPE='POLY_PWLFM'; end
if(nargin<2 || isempty(M)) M=4; end
if(nargin<3 || isempty(A)) A=[-0.001, 0, 0.001, 0.005]; end
if(nargin<4 || isempty(L)) L=ones(M,1).*(1/M); end
if(nargin<5 || isempty(K)) K=4; end
if(nargin<6 || isempty(E)) 
    E=ones(K,1);
    for n=1:K
        E(n)=E(n)./(n+2)^2;
    end
end
if(nargin<7 || isempty(p)) p=zeros(M+1,1); end
if(nargin<8 || isempty(N)) N=256; end
if(nargin<9 || isempty(epoch_length)) epoch_length=30; end % in seconds
if(nargin<10 || isempty(TV_AMP)) TV_AMP=1; end
if(nargin<11 || isempty(E_turning_points))  
    E_turning_points{1}=[5 15 25];
    E_turning_points{2}=[6 10 21];   
    E_turning_points{3}=[2 11 27];       
    E_turning_points{4}=[5 17 24];           
    E_turning_points{5}=[3 15 18];               
end
if(nargin<12 || isempty(E_tv_points))  
    E_tv_points{1}=[0.8 1.1 0.9]; 
    E_tv_points{2}=[1.1 0.8 1];     
    E_tv_points{3}=[0.7 1.1 1.2];     
    E_tv_points{4}=[1.1 1 0.7];     
    E_tv_points{5}=[0.9 0.8 1];     
end
if(nargin<13 || isempty(f0)) f0=0.12; end


y=[]; iflaw=[];


%---------------------------------------------------------------------
% SWITCHES
%---------------------------------------------------------------------
DB=0;
DBplot=0;
DBtest=0;

% Actually, don't want to create new figures:
% $$$ if(DBplot)
% $$$     openfigs=findobj('Type','figure');
% $$$     if(isempty(openfigs))
% $$$         fignum=1;
% $$$     else
% $$$         fignum=max(openfigs)+1;
% $$$     end
% $$$ end    

%---------------------------------------------------------------------
% Some defintions
%---------------------------------------------------------------------
START_FREQ=f0;
T=epoch_length/N;


%---------------------------------------------------------------------
% 1. Generate IF law
%---------------------------------------------------------------------
iflaw=zeros(N,1);
iends=floor(L.*N);
if(sum(iends)>N)  error('Segments too long.'); end


freq_last=START_FREQ;
ilast=0;
for m=1:M
    n=1:iends(m); n=n+ilast;
    iflaw(n)=A(m).*(1:length(n)).*T+freq_last;    
    
    ilast=n(end);
    freq_last=iflaw(ilast);
    
    if( strcmp(TYPE,'POLY_PWLFM') )
        iturning_points(m)=n(1);
    end
end


if( strcmp(TYPE,'POLY_PWLFM') )
    iturning_points(M+1)=n(end);
    
    iflaw=spline(iturning_points,iflaw(iturning_points),1:N);    
end



% check the limits on the IF law:
if(max(abs(iflaw))>0.5 || min(iflaw)<0)
    disp('WARNING: IF law goes out of bounds (>0.5, or <0)');
    
    if(min(iflaw)<0)
        disp('fixing....');        
        iflaw=iflaw+abs(min(iflaw))+0.01;
    elseif(max(abs(iflaw))>0.5)
        disp('fixing....');        
        iflaw=iflaw-abs(max(iflaw))-0.01;
    else 
        y=[]; iflaw=[];
        return;
    end
end



energy_fund=zeros(1,1); energy_harmonics=zeros(K,1);
energy_ratio_fund_harmonics=zeros(1,1);
%---------------------------------------------------------------------
% 2. Put IF law into signal
%---------------------------------------------------------------------
% 2a. signal:
y=cos(2.0*pi*cumtrapz(iflaw)+p(1));


if(DBplot) 
    title('IF laws');
    figure(1); clf;  
    n=0:N-1; t=n.*T;
    line_colors=get(gca,'ColorOrder');

    iflaw_estimate=get_IF(y); 
    subplot(3,1,1);  hold on;
    h=plot(t,iflaw,'-');
    set(h,'color',line_colors(1,:));
    lgd_str{1}='FUN1';
end

%---------------------------------------------------------------------
% IF time-varying amplitude
%---------------------------------------------------------------------
if(TV_AMP)
    if(isempty(E_turning_points) || isempty(E_turning_points))
        error('Need to enter a sensible value for TV-amplitude turning points');
    end
    
    tp=E_turning_points{1}(:);
    tp_amp=E_tv_points{1}(:);
    
    tv_amp=spline(tp.*(1/T),tp_amp,1:N);
    y=y(:).*tv_amp(:);
    
end
if(DBplot)
    figure(5); clf;
    title('Time-varying amplitudes');
    subplot((K+1),1,1);
    plot(t,y,'color',line_colors(1,:));
    hold on;
        
% $$$         subplot(212);
% $$$         hold on;
% $$$         plot(t,tv_amp,line_colors(1,:));
% $$$         hold on;
end


energy_fund = nanmean(abs(y) .^ 2);
if(DBtest) if_laws_test=iflaw; end



% 2b harmonics:
for kk=1:K
    CHOP_OUT=0;
    
    iflaw_harm=(kk+1).*iflaw;    
    if(max(iflaw_harm)<0.5)
        y_harmonic=sqrt(E(kk)).*cos(2.0*pi*(cumtrapz(iflaw_harm))+p(kk+1));
        
        % If time-varying amplitude:
        if(TV_AMP)
            tv_amp=spline(E_turning_points{kk+1}*(1/T),E_tv_points{kk+1},1:N);
            y_harmonic=y_harmonic(:).*tv_amp(:);
 
            if(DBplot)
                subplot(212);     hold on;
                plot(t,tv_amp,'color',line_colors(kk+1,:));
            end
        end
        if(DBplot)
            figure(5);
            subplot((K+1),1,kk+1);
            hold on;
            plot(t,y_harmonic,'color',line_colors(kk+1,:)); 
                
        end
        

        energy_harmonics(kk)=getEnorm(y_harmonic);
        
        
        y=y+y_harmonic;
    else
        disp('WARNING - Chopped out harmonic component.');        
        CHOP_OUT=1;
    end

    
    if(DBtest) if_laws_test=if_laws_test+iflaw_harm; end
    
    if(DBplot && ~CHOP_OUT)
        lgd_str{kk+1}=strcat('H:',num2str(kk));
        figure(1);        
        subplot(3,1,1);
        plot(t,iflaw_harm,'color',line_colors(kk+1,:));
        axis([t(1) t(end) 0 0.5]);        
    end
end

if(K>0)
    energy_ratio_fund_harmonics=sum(energy_harmonics)/energy_fund;
end


if(DB)
    dispVars(energy_fund,energy_harmonics,energy_ratio_fund_harmonics);
end



if(DBplot) 
    figure(1);            
    legend(lgd_str); 
    hold off;
end



if(DBplot)    
    figure(1);            
    subplot(3,1,2);
    iflaw_estimate=get_IF(get_analytic(y));
    iflaw_estimate=iflaw_estimate(1:N);     
    plot(t,iflaw,'-',t,iflaw_estimate,'-');
    axis([t(1) t(end) 0 0.5]);

    subplot(3,1,3);
    plot(t,real(y));
end

% JUST AN IDEA!!!
if(DBtest)
    iflaw_estimate=get_IF(y);     
    % aa) low pass filter:
    iflaw_estimate=filt_detrend_and_BP(iflaw_estimate,1/T,1,0,0.06,101);

% $$$     iflaw_estimate=iflaw_estimate-mean(iflaw_estimate);
    
    
    w=getTFD(iflaw_estimate,[],'sep');
    figure(311); clf;
    vtfd(w,iflaw_estimate);
    
    
% $$$     wper=getTFD(if_laws_test,[],'sep');
% $$$     figure(312); clf;
% $$$     vtfd(wper,if_laws_test);
end



if(DBplot)
    w=dtfd_sep1(y,{51,'hamm',0,1},{271,'hann'},512,512);
    figure(2); clf;
    vtfd(w,y);
end

