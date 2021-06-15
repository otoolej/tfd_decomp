%-------------------------------------------------------------------------------
% nlfm_test_signals: generate signals with nonlinear-FM laws
%
% Syntax: x = nlfm_test_signals()
%
% Inputs: 
%      - 
%
% Outputs: 
%     x - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 06-09-2019
%
% last update: Time-stamp: <2019-09-06 17:10:02 (otoolej)>
%-------------------------------------------------------------------------------
function x = nlfm_test_signals(DBplot)
if(nargin < 1 || isempty(DBplot)), DBplot = false; end


%---------------------------------------------------------------------
% signals with nonlinear FM laws:
%---------------------------------------------------------------------
M = 4;
A = [-0.001, 0, 0.001, 0.005]; 
L = ones(M,1) .* (1/M); 
num_harmonics = 0;
tv_amp = 0;
N = 256; Fs = 120;

[y1, iflaw] = gen_signal('POLY_PWLFM', M, A, L, num_harmonics, [], [], N, Fs, tv_amp);

M = 4;
A = [0.002, 0, -0.001, 0.001]; 
L = [0.3, 0.1, 0.25, 0.35];
num_harmonics = 0;
tv_amp = 0;
f0 = 0.35;

[y2, iflaw] = gen_signal('POLY_PWLFM', M, A, L, num_harmonics, [], [], N, Fs, tv_amp, ...
                         [], [], f0);


%---------------------------------------------------------------------
% and linear FM signals:
%---------------------------------------------------------------------
lfm1 = gen_LFM(N, 0.4, 0.3);
lfm2 = gen_LFM(N, 0.2, 0.4);
lfm3 = gen_LFM(N, 0.2, 0.3);
lfm3(1:50) = 0; lfm3(180:end) = 0; 

%---------------------------------------------------------------------
% combine
%---------------------------------------------------------------------
x1 = (real(y1(:)) + real(lfm1(:)))';
x2 = (real(y2(:)) + real(y1(:)))';
x3 = (real(lfm1(:)) + real(lfm2(:)))';
x4 = (real(y1(:)) + real(y2(:)) + real(lfm3(:)))';
x5 = (0.5 .* real(y1(:)) + 1 .* real(y2(:)) + 0.1 .* real(lfm3(:)))';


%---------------------------------------------------------------------
% try and replicate example in ensemble-EMD paper
%---------------------------------------------------------------------
s6 = gen_LFM(N, 0.018, 0.018, -pi/2); 
f6 = gen_LFM(N, 0.45, 0.45);
zr = zeros(N, 1);
zr(65:75) = 1; zr(120:130) = 1; zr(175:185) = 1;
f6 = (f6 .* zr.') .* 0.5;
x6 = s6 + f6;


x = [x1; x2; x3; x4; x5; x6];


if(DBplot)
    set_figure(1); 
    plot_components_(x);
end

