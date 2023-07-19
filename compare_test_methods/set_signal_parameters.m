%-------------------------------------------------------------------------------
% set_signal_parameters: default parameters for signals 
%
% Syntax: [x, x_comps, all_params] = set_signal_parameters(signal_name)
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
% last update: Time-stamp: <2023-07-18 18:19:14 (otoolej)>
%-------------------------------------------------------------------------------
function [x, x_components, Fs, all_params] = set_signal_parameters(signal_name, db_plot)
if(nargin < 2 || isempty(db_plot)), db_plot = false; end



%---------------------------------------------------------------------
% 1. load signals
%---------------------------------------------------------------------
[x_nlfm, x_nlfm_comps, x_lfm_comps] = nlfm_test_signals(false);
Fs = 30;
N = 256;

% find length of signal (for setting default parameters)
if(strcmp(signal_name, 'bat'))
    d = load('batsignal');
    N = length(d.batsignal);

elseif(strcmp(signal_name(1:3), 'nnl'))
    idx = str2num(signal_name(end));
    N = length(x_nlfm{idx});
    
elseif(strcmp(signal_name, 'noise'))
    N = 128;
end


%---------------------------------------------------------------------
% 2. default parameters for all methods
%---------------------------------------------------------------------
[xtfd, tvfilt, ssst, wsst, vmd, vncmd, tvemd, ncme, msst] = set_default_params_all_methods(N);



switch signal_name
    
  case {'2tone1'}
    [x, x_components] = get_2tone_test_signal(false, 25, 40); 
    
    N = 1024;
    vmd.alpha = 1;
    vmd.tau = 0;
    vmd.DC = 0;
    vmd.init = 1;
    vmd.tol = 1e-9;

    vncmd.estIF = [307 .* ones(1, N); 410 .* ones(1, N)];
    vncmd.alpha = 1;
    vncmd.beta = 8e-6;
    vncmd.tol = 1e-9;
    ssst.window = chebwin(67, 100);

    tvemd.bwr = 0.05;
    
    tvfilt = decomp_params(N, 'tvfilt');
    xtfd = decomp_params(N, 'xtfd');

    l_lag = make_odd(ceil(sqrt(N) * 2));
    l_dopp = make_odd(ceil(sqrt(N)));
    tvfilt = tvfilt.set_lag_kernel(l_lag, {'dolph', 100});
    tvfilt = tvfilt.set_dopp_kernel(l_dopp, {'hamm'});
    tvfilt.qtfd_max_thres = [];

    l_dopp = make_odd(N / 2);
    xtfd = xtfd.set_lag_kernel(l_lag, {'dolph', 100});
    xtfd = xtfd.set_dopp_kernel(l_dopp, {'hamm'});
    xtfd.Nfreq = N * 8;
    
    
  case {'2tone2'}    
    [x, x_components] = get_2tone_test_signal(true, 25, 10); 
    
    N = 1024;    
    vmd.alpha = 1;
    vmd.tau = 0;
    vmd.DC = 0;
    vmd.init = 1;
    vmd.tol = 1e-9;

    vncmd.estIF = [82 .* ones(1, N); 409 .* ones(1, N)];
    vncmd.alpha = 1;
    vncmd.beta = 8e-6;
    vncmd.tol = 1e-9;
    
    tvemd.bwr = 0.05;
    
    ssst.window = chebwin(33, 100);
    
    tvfilt = decomp_params(N, 'tvfilt');
    xtfd = decomp_params(N, 'xtfd');

    l_lag = make_odd(ceil(sqrt(N) * 2));
    l_dopp = make_odd(ceil(sqrt(N)));
    tvfilt = tvfilt.set_lag_kernel(l_lag, {'dolph', 100});
    tvfilt = tvfilt.set_dopp_kernel(l_dopp, {'hamm'});
    tvfilt.qtfd_max_thres = [];

    l_dopp = make_odd(N / 2);
    xtfd = xtfd.set_lag_kernel(l_lag, {'dolph', 100});
    xtfd = xtfd.set_dopp_kernel(l_dopp, {'hamm'});
    xtfd.Nfreq = N * 8;
    

  case {'nnlfm4'}
    % same as previous but with noise
    idx = str2num(signal_name(end));    
    x = x_nlfm{idx};
    x_components = x_nlfm_comps{idx};
    % N_components = 50;

    d = load('data/test_signals/ffgn_1_02_1_512_0_signal.mat');
    n = d.x(1:256);
    n = n .* (0.6081);
    
    % rng(1);
    % n = randn(size(x)) * 0.55;
    fprintf('SNR=%g (dB)\n', 10 * log10(sum(abs(x) .^ 2) / sum(abs(n) .^ 2)));
    x = x + n;

    % xtfd = xtfd.set_dopp_kernel(67, {'hamm'});
    % xtfd = xtfd.set_lag_kernel( 101, {'dolph', 100});
    % default is fine but to harmonize with 'tvfilt' then:
    % xtfd.min_if_length = 64;

    l_lag = make_odd(ceil(sqrt(N) * 6));    
    l_dopp = make_odd(ceil(sqrt(N) * 4));
    
    tvfilt = tvfilt.set_dopp_kernel(l_dopp, {'hamm'});
    tvfilt = tvfilt.set_lag_kernel(l_lag, {'dolph', 100});
    xtfd = xtfd.set_dopp_kernel(l_dopp, {'hamm'});
    xtfd = xtfd.set_lag_kernel(l_lag, {'dolph', 100});

    xtfd.Nfreq = N;    
    % tvfilt.qtfd_max_thres = [];
    % tvfilt.min_if_length = 64;
    % disp(tvfilt)
    % x_components{4} = n;
    % N_components = 4;


    vncmd.estIF = [80 * ones(1, length(x)); 100 * ones(1, length(x)); ...
                   150 * ones(1, length(x))];
    vncmd.alpha = 1e-2;
    
    
  case 'noise'
    d = load('data/test_signals/ffgn_1_02_1_512_0_signal.mat');
    % x = d.x(1:256);
    x = d.x(1:128);    
    
    N_components = 20;
    for n = 1:N_components
        x_components{n} = x;
    end
    
    xtfd.min_if_length = 16;
    tvfilt.min_if_length = 16;    
    xtfd.Nfreq = 128;
    
    N = length(x);
    l_lag = make_odd(ceil(N / 4));    
    l_dopp = make_odd(ceil(N / 2));
    xtfd = xtfd.set_dopp_kernel(l_dopp, {'hamm'});
    xtfd = xtfd.set_lag_kernel(l_lag, {'dolph', 100});
    tvfilt = tvfilt.set_dopp_kernel(l_dopp, {'hamm'});
    tvfilt = tvfilt.set_lag_kernel(l_lag, {'dolph', 100});

    tvemd.bwr = 0.1;
    
    ssst.window = chebwin(63, 100);
    vncmd.estIF = [50 * ones(1, length(x)); 80 * ones(1, length(x)); ...
                   100 * ones(1, length(x)); 200 * ones(1, length(x))];
    vncmd.alpha = 0.01;
    vncmd.beta = 1;
    tvemd.bwr = 0.1;
    vmd.alpha = 1;
    msst.hlength = 60;
    ncme.estIF = [10 * ones(1, length(x)); 20 * ones(1, length(x)); ...
                  50 * ones(1, length(x))];
    ncme.beta = 1;
    ncme.lambda = 0.01;
    



  case {'bat'}
    N_components = 50;
    d = load('batsignal');
    x = d.batsignal;
    % dummy components:
    for p = 1:N_components
        x_components{p} = x;
    end

    N = length(x);
    l_lag = make_odd(ceil(sqrt(N) * 2));
    l_dopp = make_odd(ceil(sqrt(N) * 4));

    xtfd = xtfd.set_lag_kernel(l_lag, {'dolph', 100});
    xtfd = xtfd.set_dopp_kernel(l_dopp, {'hamm'});    
    tvfilt = tvfilt.set_lag_kernel(l_lag, {'dolph', 100});
    tvfilt = tvfilt.set_dopp_kernel(l_dopp, {'hamm'});    
    xtfd.Nfreq = N;        

    ssst.window = chebwin(63, 100);
    vncmd.estIF = [50 * ones(1, length(x)); 80 * ones(1, length(x)); ...
                   100 * ones(1, length(x)); 200 * ones(1, length(x))];
    vncmd.alpha = 0.001;
    vncmd.beta = 1;
    tvemd.bwr = 0.1;
    vmd.alpha = 1;
    msst.hlength = 60;

    ncme.estIF = [50 * ones(1, length(x)); 30 * ones(1, length(x)); ...
                  60 * ones(1, length(x)); 20 * ones(1, length(x))];
    ncme.lambda = 1;
    ncme.beta = 1e-6;


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



function [xtfd, tvfilt, ssst, wsst, vmd, vncmd, tvemd, ncme, msst] = ...
    set_default_params_all_methods(N)
%---------------------------------------------------------------------
% default parameters for all methods
%---------------------------------------------------------------------
if(nargin < 1 || isempty(N)), N = 256; end


%---------------------------------------------------------------------
% 2. default parameters for all methods
%---------------------------------------------------------------------
% xTFD method
xtfd = decomp_params(N, 'xtfd');


% TV filtering method
tvfilt = decomp_params(N, 'tvfilt');


% STFT synchrosqueeze:
ssst.window = chebwin(63, 100);
% wavelet synchrosqueeze
wsst.window = 'bump';

% VMD
vmd.alpha = 1e-4;
vmd.tau = 0;
vmd.DC = 0;
vmd.init = 1;
vmd.tol = 1e-9;

% VNCMD
vncmd.estIF = 1;
vncmd.alpha = 1e-6;
vncmd.beta = 5e-6;
vncmd.var = 0;
vncmd.tol = 1e-9;

% TV-EMD
tvemd.bwr = 0.01;


% VNCMD
ncme.estIF = floor(N / 2) .* ones(1, N);
ncme.lambda = 100;
ncme.beta = 5e-6;
ncme.tol = 1e-9;


% multisychronsqueezing transform
msst.hlength = 50;
