%-------------------------------------------------------------------------------
% select_signal_withparams: 
%
% Syntax: [x, x_comps, all_params] = select_signal_withparams(signal_name)
%
% Inputs: 
%     signal_name - 
%
% Outputs: 
%     [x, x_comps, all_params] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 24-09-2021
%
% last update: Time-stamp: <2021-12-06 11:01:02 (otoolej)>
%-------------------------------------------------------------------------------
function [x, x_components, Fs, all_params] = select_signal_withparams(signal_name, db_plot)
if(nargin < 2 || isempty(db_plot)), db_plot = false; end



%---------------------------------------------------------------------
% 1. load signals
%---------------------------------------------------------------------
[x_nlfm, x_nlfm_comps, x_lfm_comps] = nlfm_test_signals(false);
Fs = 30;
N = 256;


%---------------------------------------------------------------------
% 2. default parameters for all methods
%---------------------------------------------------------------------
[xtfd, tvfilt, ssst, wsst, vmd, vncmd, tvemd, ncme, msst] = set_default_params_all_methods(N);



switch signal_name
    %---------------------------------------------------------------------
    % linear FM signals
    %---------------------------------------------------------------------

  case {'lfm1', 'lfm2', 'lfm3'}
    idx = str2num(signal_name(end));
    x = x_lfm_comps{idx};
    x_components = {x_lfm_comps{idx}};
    N_components = 1;

    vmd.alpha = 1e-4;
    xtfd.lag_kernel = {67, 'dolph', 100};
    ssst.window = chebwin(61, 100);
    wsst.window = 'amor';
    vncmd.estIF = 100 * ones(1, length(x));
    vncmd.alpha = 8;
    vncmd.beta = 1e-6;
    tvfilt.wx = 67;

    ncme.estIF = vncmd.estIF;
    ncme.lambda = 10;
    ncme.beta = 3e-6;
    ncme.tol = 1e-9;


    % xtfd.lag_kernel = {11, 'dolph', 100};
    % xtfd.doppler_kernel = @(x){27, 'hamm', 100};        
    

    switch signal_name
      case 'lfm3'
        vncmd.alpha = 1;
        vncmd.beta = 1;
    end

    % N_components = 2;

    
  case {'lfm4', 'lfm5', 'lfm6'}
    idx = str2num(signal_name(end));
    % x_lfm_comps{idx}{2} = x_lfm_comps{idx}{2} ./ 8;
    x = sum(cat(1, x_lfm_comps{idx}{:}), 1);
    x_components = x_lfm_comps{idx};
    N_components = 2;

    % x_components = x_lfm_comps{idx}(1);
    % x = x_components{1};
    % N_components = 1;

    if(idx == 4)
        xtfd.lag_kernel = {63, 'dolph', 100};        
        ssst.window = chebwin(63, 100);

        vmd.alpha = 10;
        % sensitive to this parameter
        vncmd.estIF = [50 * ones(1, length(x)); 100 * ones(1, length(x))];
        vncmd.alpha = 1;
        vncmd.beta = 1;

        ncme.estIF = vncmd.estIF;
        ncme.beta = 1;
        
    elseif(idx == 5)
        xtfd.lag_kernel = {127, 'dolph', 100};
        xtfd.doppler_kernel = @(x){27, 'hamm', 100};        
        ssst.window = chebwin(227, 100);
        tvfilt.wx = 67;

        vmd.alpha = 100;
        % sensitive to this parameter
        vncmd.estIF = [92 * ones(1, length(x)); 113 * ones(1, length(x))];
        vncmd.alpha = 5;
        vncmd.beta = 1;

        ncme.estIF = vncmd.estIF;
        ncme.beta = 1;

    elseif(idx == 6)
        % xtfd.lag_kernel = {67, 'dolph', 100};
        % tvfilt.wx = 11;


        vmd.alpha = 100;
        % sensitive to this parameter
        vncmd.estIF = [100 * ones(1, length(x)); 100 * ones(1, length(x))];
        vncmd.alpha = 5;
        vncmd.beta = 1e-6;

        ncme.estIF = vncmd.estIF;
        ncme.beta = 1e-6;
        
    end

  case {'nlfm1', 'nlfm2', 'nlfm3'}
    idx = str2num(signal_name(end));    
    x = x_nlfm{idx};
    x_components = x_nlfm_comps{idx};
    N_components = 2;

    vncmd.estIF = [100 * ones(1, length(x)); 120 * ones(1, length(x))];
    
  case 'msst2'
    % same as below but without noise
    SampFreq = 100;
    t = 1/SampFreq : 1/SampFreq : 4;
    x_components = {sin(2*pi*(40*t + 1*sin(1.5*t))); sin(2*pi*(17*t + 6*sin(1.5*t)))};
    x = sum(cat(1, x_components{:}));
    N_components = 2;
    
    vncmd.estIF = [100 * ones(1, length(x)); 120 * ones(1, length(x))];
    vncmd.alpha = 1;
    vncmd.beta = 1;

    % xtfd.lag_kernel = {127, 'dolph', 100};
    xtfd.doppler_kernel = @(x){27, 'hamm', 100};        
    ssst.window = chebwin(227, 100);
    tvfilt.wx = 27;
    % tvfilt.M = 21;
    

  case 'msst_test'
    % from the multi-SST paper:
    d = load('Sig_noise.mat');
    x = d.Sig_noise(12, :);
    
    SampFreq = 100;
    t = 1/SampFreq : 1/SampFreq : 4;
    x_components = {sin(2*pi*(40*t + 1*sin(1.5*t))); sin(2*pi*(17*t + 6*sin(1.5*t)))};
    N_components = 2;

    xtfd.lag_kernel = {63, 'dolph', 100};        
    ssst.window = chebwin(63, 100);

    vmd.alpha = 10;
    % sensitive to this parameter
    vncmd.estIF = [50 * ones(1, length(x)); 100 * ones(1, length(x))];
    vncmd.alpha = 1;
    vncmd.beta = 1;
    
    

  case {'vmd_test'}
    % from the VMD paper:
    T = 1000;
    fs = 1/T;
    t = (1:T)/T;
    % center frequencies of components
    f_1 = 2;
    f_2 = 24;
    f_3 = 288;

    % modes
    x_components{1} = (cos(2*pi*f_1*t));
    x_components{2} = 1/4*(cos(2*pi*f_2*t));
    x_components{3} = 1/16*(cos(2*pi*f_3*t));
    N_components = 3;
    x = sum(cat(1, x_components{:})) + 0.1 * randn(size(x_components{1}));
    Fs = 100;
    

    xtfd.lag_kernel = {127, 'dolph', 100};        
    

  case {'nlfm4', 'nlfm5', 'nlfm6'}
    idx = str2num(signal_name(end));    
    x = x_nlfm{idx};
    x_components = x_nlfm_comps{idx};
    N_components = 3;

    vncmd.estIF = [80 * ones(1, length(x)); 100 * ones(1, length(x)); ...
                  150 * ones(1, length(x))];

  case {'nlfm7'}
    % from the VNCMD paper
    idx = str2num(signal_name(end));        
    x = x_nlfm{idx}{:};
    x_components = x_nlfm_comps{idx};
    N_components = 1;
    Fs = 200;

    x = x +  (0.8 .* randn(1, length(x)));
    vncmd.estIF = [100 * ones(1, length(x))];

    xtfd.lag_kernel = {223, 'dolph', 100};
    ssst.window = chebwin(223, 100);
    
    
  case {'nlfm8'}
    % from the VNCMD paper
    idx = str2num(signal_name(end));        
    x = x_nlfm{idx};
    x_components = x_nlfm_comps{idx};
    N_components = 2;
    Fs = 200;
    
    xtfd.lag_kernel = {23, 'dolph', 100};
    % xtfd.doppler_kernel = @(x){27, 'hamm', 100};            
    ssst.window = chebwin(63, 100);
    vncmd.estIF = [50 * ones(1, length(x)); 100 * ones(1, length(x))];


  case 'noise'
    d = load('data/ffgn_1_02_1_512_0_signal.mat');
    x = d.x;
    % x = ffgn(1, 0.2, 1, 512, 0);
    % detrend:
    % y = filt_butterworth(x, 1, 0.05, [], 7, false);
    % x = x - y;
    x_components = {x, x, x, x, x, x, x, x, x, x, ...
                    x, x, x, x, x, x, x, x, x, x};
    N_components = 20;
    xtfd.lag_kernel = {17, 'dolph', 100};
    % xtfd.doppler_kernel = @(x){33, 'hamm'};
    xtfd.delta_search_freq = 1000;
    xtfd.max_no_peaks = 200;
    % xtfd.Nfreq = 512;
    tvfilt.wx = 17;

    

  case {'bat'}
    N_components = 15;
    d = load('batsignal');
    x = d.batsignal;
    % dummy components:
    for p = 1:N_components
        x_components{p} = x;
    end

    xtfd.lag_kernel = {63, 'dolph', 100};
    xtfd.doppler_kernel = @(x){63, 'hamm'};    
    ssst.window = chebwin(63, 100);
    vncmd.estIF = [50 * ones(1, length(x)); 80 * ones(1, length(x)); ...
                   100 * ones(1, length(x)); 200 * ones(1, length(x))];
    vncmd.alpha = 0.001;
    vncmd.beta = 1;
    tvemd.bwr = 0.1;
    vmd.alpha = 1;
    msst.hlength = 60;
    tvfilt.wx = 31;

  otherwise
    error('what signal?');

end
x = x(:);


all_params(1) = struct('method', 'xtfd', 'params', xtfd);
all_params(2) = struct('method', 'tvfilt', 'params', tvfilt);
all_params(3) = struct('method', 'ssst', 'params', ssst);
all_params(4) = struct('method', 'wsst', 'params', wsst);
all_params(5) = struct('method', 'vmd', 'params', vmd);
all_params(6) = struct('method', 'vncmd', 'params', vncmd);
all_params(7) = struct('method', 'tvemd', 'params', tvemd);
all_params(8) = struct('method', 'ixtfd', 'params', xtfd);
all_params(9) = struct('method', 'msst', 'params', msst);
all_params(10) = struct('method', 'ncme', 'params', ncme);



if(db_plot)
    plot_components(x, x_components, 1, 222, length(x_components), false);
    title('test signals');
end
