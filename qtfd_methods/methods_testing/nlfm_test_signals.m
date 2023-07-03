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
% last update: Time-stamp: <2021-10-01 18:35:31 (otoolej)>
%-------------------------------------------------------------------------------
function [x, x_components, x_lfm_components] = nlfm_test_signals(db_plot)
if(nargin < 1 || isempty(db_plot)), db_plot = false; end


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
x_lfm = zeros(3, N);
x_lfm(1, :) = gen_LFM(N, 0.3, 0.3);
x_lfm(2, :) = gen_LFM(N, 0.2, 0.4);
x_lfm(3, :) = gen_LFM(N, 0.2, 0.3);
x_lfm(3, 1:50) = 0; 
x_lfm(3, 180:end) = 0;

x_lfm_components = num2cell(x_lfm.', 1);
x_lfm_components{4} = {gen_LFM(N, 0.1, 0.15), 0.25 .* x_lfm(3, :)};
x_lfm_components{5} = {gen_LFM(N, 0.18, 0.21), 2.75 .* gen_LFM(N, 0.22, 0.25)};

x_stat = gen_LFM(N, 0.2, 0.2);

mask1 = zeros(1, N);
mask1(1:floor(N / 3)) = 1;
x11 = x_stat .* mask1;
x22 = x_stat .* fliplr(mask1);
x_lfm_components{6} = {x11, x22};

%---------------------------------------------------------------------
% combine
%---------------------------------------------------------------------
x_components{1} = {0.2 .* real(y1), real(x_lfm(1, :))};
x_components{2} = {real(y2), real(y1)};
x_components{3} = {real(x_lfm(1, :)), real(x_lfm(2, :))};
x_components{4} = {real(y1), real(y2), real(x_lfm(3, :))};
x_components{5} = {0.5 .* real(y1), 1 .* real(y2), 0.1 .* real(x_lfm(3, :))};

% from 'Multisynchrosqueezing Transform' (Yu et al, 2019)
SampFreq = 100;
t = 1/SampFreq : 1/SampFreq : 4;
x_components{6} = {sin(2*pi*(17*t + 6*sin(1.5*t))), sin(2*pi*(40*t + 1*sin(1.5*t)))};



%---------------------------------------------------------------------
% from VNCMD paper
%---------------------------------------------------------------------
fs = 2000;
t = 1 / fs:1 / fs:1;
x_components{7} = {cos(2*pi*(1.5 + 150*t - 100*(t.^2) + 416*(t.^3) - 200*(t.^4)))};
x_components{8} = {cos(390 * pi * t + 0.2 * sin(10 * pi * t)), ...
                   cos(400 * pi * t + 0.2 * sin(10 * pi * t))};
% which is equal to this:
x_test = 2 * cos(-5 * pi .* t) .* cos((790 / 2) * pi .* t + 0.2 * sin(10 * pi .* t));


for n = 1:length(x_components)
    if(length(x_components{n}) > 1)
        x{n} = sum(cat(1, x_components{n}{:}));
    else
        x{n} = x_components{n};
    end
end


%---------------------------------------------------------------------
% try and replicate example in ensemble-EMD paper
%---------------------------------------------------------------------
s6 = gen_LFM(N, 0.018, 0.018, -pi/2); 
f6 = gen_LFM(N, 0.45, 0.45);
zr = zeros(N, 1);
zr(65:75) = 1; zr(120:130) = 1; zr(175:185) = 1;
f6 = (f6 .* zr.') .* 0.5;
x6 = s6 + f6;
L = length(x_components);
x{L + 1} = x6;
x_components{L + 1} = {s6, f6};



if(db_plot)
    set_figure(1);
    plot_components_(cell2mat(x'));
end

