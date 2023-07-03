%-------------------------------------------------------------------------------
% plot_2tone_test_signals: 
%
% Syntax:  = plot_2tone_test_signals()
%
% Inputs: 
%      - 
%
% Outputs: 
%      - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 11-02-2022
%
% last update: Time-stamp: <2022-04-24 14:17:30 (otoolej)>
%-------------------------------------------------------------------------------
function plot_2tone_test_signals(print_)
if(nargin < 1 || isempty(print_)), print_ = false; end



%---------------------------------------------------------------------
% parameters for the signals
%---------------------------------------------------------------------
L_grid = 50;
L_freq = L_grid;
L_amp = L_grid;
    
N = 1024;
f_max = 0.8;
f_start = 0.01;

k0 = linspace(f_start, f_max - f_start, L_freq);
amp_ratio = 10 .^ linspace(-2, 2, L_amp);
Fs = 32;
T = 1 / Fs;
t = T:T:(N * T);

% frequencies and ratio:
f0 = (f_max + f_start) * (Fs / 2);
f1 = k0 .* (Fs / 2);
f_ratio = f1 ./ f0;

% first component:
x_comp{1} = cos(2 * pi * f0 .* t);

k = 12;
m = 30;
        
x_comp{2} = amp_ratio(m) .* cos(2 * pi * f1(k) .* t);


postproc_plot_components(sum(cat(1,x_comp{:})), x_comp, 900, false);
pp = get(gcf, 'position');
set(gcf, 'position', [pp(1:2) 600  240]);
% set(gcf, 'position', [pp(1:2) 770  140]);

if(print_)
    fname = ['pics/twotone_test1_eg1_v3'];
    % print([fname '.svg'], '-dsvg');
    print2eps([fname '.eps']);
end


k = 5;
m = 40;

x_comp{2} = amp_ratio(m) .* cos(2 * pi * f1(k) .* t);
time_mask = true;
if(time_mask)
    mask1 = zeros(1, N);
    mask1(1:floor(N / 3)) = tukeywin(floor(N/3), 0.2);
    x_comp{1} = x_comp{1} .* mask1;
    x_comp{2} = x_comp{2} .* fliplr(mask1);    
end

postproc_plot_components(sum(cat(1,x_comp{:})), x_comp, 800, false);
pp = get(gcf, 'position');
set(gcf, 'position', [pp(1:2) 570  240]);
% set(gcf, 'position', [pp(1:2) 770  140]);

if(print_)
    fname = ['pics/twotone_test2_eg1_v3'];
    % print([fname '.svg'], '-dsvg');
    print2eps([fname '.eps']);
end
