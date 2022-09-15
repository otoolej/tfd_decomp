%-------------------------------------------------------------------------------
% set_default_params_all_methods: 
%
% Syntax: [xtfd, tvfilt, ssst, wsst, vmd, vncmd, tvemd, efd] = set_default_params_all_methods()
%
% Inputs: 
%      - 
%
% Outputs: 
%     [xtfd, tvfilt, ssst, wsst, vmd, vncmd, tvemd, efd] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 24-09-2021
%
% last update: Time-stamp: <2021-11-30 17:54:35 (otoolej)>
%-------------------------------------------------------------------------------
function [xtfd, tvfilt, ssst, wsst, vmd, vncmd, tvemd, ncme, msst] = set_default_params_all_methods(N)
if(nargin < 1 || isempty(N)), N = 256; end



%---------------------------------------------------------------------
% 2. default parameters for all methods
%---------------------------------------------------------------------
% xTFD method
xtfd = decomp_params;
xtfd.lag_kernel = {63, 'dolph', 100};

% TV filtering method
tvfilt = tvfilt_params(N);


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
