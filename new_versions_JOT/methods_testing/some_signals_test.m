%-------------------------------------------------------------------------------
% some_signals_test: testing the TF decomposition method with some test signals
%
% Syntax: []=some_signals_test(TEST_SIGNAL)
%
% Inputs: 
%     TEST_SIGNAL - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University of Deusto
% Started: 24-01-2013
%-------------------------------------------------------------------------------
function []=some_signals_test(TEST_SIGNAL, PRINT_)
if(nargin<1 || isempty(TEST_SIGNAL)) TEST_SIGNAL=3; end
if(nargin<2 || isempty(PRINT_)) PRINT_=0; end



% load the signals:
[x, x_lfm] = nlfm_test_signals(false);
% .. and some parameters
% decomp_parameters;

DBplot = 1;

Fs=30;

%---------------------------------------------------------------------
% set the TFD kernel parameters here, as same length signals for all
%---------------------------------------------------------------------
N = size(x, 2);
tfd_params = {{floor(N/2)-1,'gauss',1,1},{floor(N/2),'hann',0,0,'n'}};


switch TEST_SIGNAL
  case -1
    % 1. same as signal '0' (LFM only) but this time estimate without phase value    
    x_test = x_lfm(2, :);
    N_remove=1; 
    
    [x_clean, x_remove, s_est, ip_tracks] = extract_components(x_test, Fs, N_remove, ...
                                                      tfd_params, DBplot, 0, [], false);
    
  case 0
    % 1. LFM only
    x_test = x_lfm(2, :);
    N_remove=1; 
    
    % x_test = linspace(1, 3, length(x_test)) .* x_test;
    % x_test(1:40) = 0;
    % % x_test(100:end) = 0;    
    % x_test(200:end) = 0;

    [x_clean, x_remove, s_est] = extract_components(x_test, Fs, N_remove, tfd_params, 1, PRINT_);

  case 1
    % 1. LFM with nonlinear-FM signal
    N_remove=2; 
    x_test = x(1, :);

    % plot_components(x1,x1_separate,Fs,17,2,0);
    % x1 = x1(1:64); % zeros(1, 256)];

    % N = length(x_test);    
    % tfd_params={{floor(N/2)-1,'gauss',5,1},{floor(N/2),'hann',0,0,'n'}};
% $$$       tfd_params=[];
    
    % x_test(1:89) = 0;
    % x_test(169:end) = 0;    
    
    
    [x_clean,x_remove,s_est]=extract_components(x_test, Fs, N_remove, tfd_params, DBplot);
    
  case 2
    % 2. 2 nonlinear-FM signal
    x_test = x(2, :);
    N_remove = 2; 
    % x2_separate{2}=real(y1); x2_separate{1}=real(y2);
    % plot_components(x2,x2_separate,Fs,17,2,0);
    
    [x_clean,x_remove,s_est]=extract_components(x_test,Fs,N_remove,tfd_params,1);
    
  case 3
    % 2. 2 nonlinear-FM signal
    x_test = x(3, :);
    N_remove=3; 
    % x4_separate{2}=real(y1); x4_separate{1}=real(y2); x4_separate{3}=real(lfm3);
    % plot_components(x4,x4_separate,Fs,17,3,1);
    
    % tfd_params={{floor(N/2)-1,'gauss',5,1},{floor(N/2),'hann',0,0,'n'}};

    [x_clean,x_remove,s_est]=extract_components(x_test,Fs,N_remove,tfd_params,1);
    
      
  case 4
      % 2. 3 nonlinear-FM signal, same as '3' but with different amplitudes
      N_remove=3; 
      x_test = x(4, :);
      % x4_separate{2}=0.5.*real(y1); x4_separate{1}=1.*real(y2); x4_separate{3}=0.1.*real(lfm3);
      % plot_components(x4_la,x4_separate,Fs,17,3,1);
      

      [x_clean,x_remove,s_est]=extract_components(x_test,Fs,N_remove,tfd_params,1);
      
      % plot_components(x4_la,s_est,Fs,18,3,1);
      
  case 5
    % 2. 3 nonlinear-FM signal, same as '3' but with different amplitudes
    %    
    %    use iterative method here
    %
    N_remove=3; 
    x_test = x(4, :);
    %     x4_separate{2}=0.5.*real(y1); x4_separate{1}=real(y2); x4_separate{3}=0.1.*real(lfm3);
    % % $$$       plot_components(x4_la,x4_separate,Fs,17,3,1);
    
    %     tfd_params={{floor(N/2)-1,'gauss',5,1},{floor(N/2),'hann',0,0,'n'}};

    [x_clean,s_est]=iterate_decomposition(x_test,Fs,1,4,tfd_params);

    % plot_components(x4_la,s_est,Fs,18,3,1);

  case 6
    % 2. 2 nonlinear-FM signal
    N_remove=4; 
    x_test = x(1, :);
    % x5_separate{1}=s5; x5_separate{2}=f5;
    % plot_components(x5,x5_separate,Fs,17,2,0);
    
    [x_clean,x_remove,s_est]=extract_components(x_test,Fs,N_remove, tfd_params, 1);
end




function plot_for_print(ip_tracks)
%---------------------------------------------------------------------
% left-over from old code;
% 
% plots IF and IP functions
%---------------------------------------------------------------------
% plot phase function
set_figure(99); 
plot(ip_tracks{1}(:,1),ip_tracks{1}(:,2));
xlim([0 ip_tracks{1}(end,1)]);        
ylim([-.55 0.05]);
u=get(gcf,'Position'); set(gcf,'Position',[u(1), u(2), u(3), u(4)/2]);        
set(gca,'xtick',[]);
set(gca,'ytick',[-.5,0]);        
set(gca,'FontName',art_params.FONT_NAME,'FontSize',art_params.FONT_SIZE);
xlabel('TIME','FontName',art_params.FONT_NAME,'FontSize', ...
       art_params.FONT_SIZE);
ylabel('PHASE','FontName',art_params.FONT_NAME,'FontSize', ...
       art_params.FONT_SIZE);

% TFD of estimated component:
% set_figure(304); title([]);
% set(gca,'clim',[0,1]);
% c=colorbar;
% set(c,'ytick',[0,1]);
% set(gca,'FontName',art_params.FONT_NAME,'FontSize',art_params.FONT_SIZE);
% % turn off axis:
% set(gca,'xtick',[]);  set(gca,'ytick',[]);
% % $$$         set(gca,'ytick',[2,4,6,8]);        
% ylabel('TIME','FontName',art_params.FONT_NAME,'FontSize', ...
%        art_params.FONT_SIZE);
% xlabel('FREQUENCY','FontName',art_params.FONT_NAME,'FontSize', ...
%        art_params.FONT_SIZE); 

% % time-domain signal:
% set_figure(1);
% u=get(gcf,'Position'); set(gcf,'Position',[u(1), u(2), u(3), u(4)/2]);        
% set(gca,'xtick',[]);
% ylim([-1.1,2]);
% set(gca,'ytick',[-1,0,1]);
% set(gca,'FontName',art_params.FONT_NAME,'FontSize',art_params.FONT_SIZE);
% xlabel('TIME','FontName',art_params.FONT_NAME,'FontSize', ...
%        art_params.FONT_SIZE);
% h=legend('signal','component 1');
% set(h,'FontName',art_params.FONT_NAME,'FontSize',art_params.FONT_SIZE);        

    
% set_figure(2);
% u=get(gcf,'Position'); set(gcf,'Position',[u(1), u(2), u(3), u(4)/2]);        
% set(gca,'xtick',[]);
% ylim([-1.1,2]);
% set(gca,'ytick',[-1,0,1]);
% set(gca,'FontName',art_params.FONT_NAME,'FontSize',art_params.FONT_SIZE);
% xlabel('TIME','FontName',art_params.FONT_NAME,'FontSize', ...
%        art_params.FONT_SIZE);
% h=legend('residual signal');
% set(h,'FontName',art_params.FONT_NAME,'FontSize',art_params.FONT_SIZE);


