function [x, x_comp] = get_2tone_test_signal(time_mask, m, k, method_str)
%---------------------------------------------------------------------
% generate an example of the 2-tone signal from the '2-tone tests'
%---------------------------------------------------------------------
if(nargin < 1 || isempty(time_mask)), time_mask = false; end
if(nargin < 2 || isempty(m)), m = 20; end
if(nargin < 3 || isempty(k)), k = 20; end
if(nargin < 4 || isempty(method_str)), method_str = 'xtfd'; end



N = 1024;
L_freq = 50;
L_amp = L_freq;
f_max = 0.8;
f_start = 0.01;
if(strcmp(upper(method_str), 'EMD'))
    f_max = 0.2;
    f_start = 0.0025;
end

k0 = linspace(f_start, f_max - f_start, L_freq);
% amp_ratio = linspace(0.01, 100, L_amp);
amp_ratio = 10 .^ linspace(-2, 2, L_amp);
Fs = 32;
T = 1 / Fs;
t = T:T:(N * T);

% frequencies and ratio:
f0 = (f_max + f_start) * (Fs / 2);
f1 = k0 .* (Fs / 2);
f_ratio = f1 ./ f0;


x_comp{1} = cos(2 * pi * f0 .* t);

if(time_mask)
    mask1 = zeros(1, N);
    mask1(1:floor(N / 3)) = tukeywin(floor(N/3), 0.2);
    x_comp{1} = x_comp{1} .* mask1;
end

x_comp{2} = amp_ratio(m) .* cos(2 * pi * f1(k) .* t);
if(time_mask)
    x_comp{2} = x_comp{2} .* fliplr(mask1);
end

fprintf('|k=%d |m=%d | f_ratio=%g | k0=%g | amp_ratio=%g |\n', ...
        k, m, f_ratio(k), k0(k), amp_ratio(m));


x = x_comp{1} + x_comp{2};


set_figure(98); 
plot(x);
plot(x_comp{1} - 5);
plot(x_comp{2} - 10);


