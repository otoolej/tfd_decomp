%-------------------------------------------------------------------------------
% sep_2tones_test: 
%
% Syntax: [] = sep_2tones_test(method)
%
% Inputs: 
%     method - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 14-09-2021
%
% last update: Time-stamp: <2023-07-07 14:27:35 (otoolej)>
%-------------------------------------------------------------------------------
function re = sep_2tones_test(method, time_mask, save_, db_plot, L_grid)
if(nargin < 1 || isempty(method)), method = 'xTFD'; end
if(nargin < 2 || isempty(time_mask)), time_mask = false; end
if(nargin < 3 || isempty(save_)), save_ = false; end
if(nargin < 4 || isempty(db_plot)), db_plot = false; end
if(nargin < 5 || isempty(L_grid)), L_grid = 50; end



%---------------------------------------------------------------------
% parameters for the signals
%---------------------------------------------------------------------
L_freq = L_grid;
L_amp = L_grid;

    
N = 1024;
f_max = 0.8;
f_start = 0.01;
if(strcmp(upper(method), 'EMD'))
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


%---------------------------------------------------------------------
% parameters for the methods
%---------------------------------------------------------------------
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


if(~time_mask)
    ssst.window = chebwin(67, 100);
else
    ssst.window = chebwin(33, 100);
end
vmd.alpha = 1;
vmd.tau = 0;
vmd.DC = 0;
vmd.init = 1;
vmd.tol = 1e-9;

ncme.estIF = [250 .* ones(1, N); 400 .* ones(1, N)];
ncme.lambda = 1000;
ncme.beta = 8e-6;
ncme.tol = 1e-9;

tvemd.bwr = 0.05;

all_params(1) = struct('method', 'xtfd', 'params', xtfd);
all_params(2) = struct('method', 'tvfilt', 'params', tvfilt);
all_params(3) = struct('method', 'tvemd', 'params', tvemd);
all_params(4) = struct('method', 'vmd', 'params', vmd);
all_params(5) = struct('method', 'ssst', 'params', ssst);
all_params(6) = struct('method', 'ncme', 'params', ncme);



p_st = all_params(find(strcmp({all_params.method}, lower(method))));
if(isempty(p_st))
    params = [];
else
    params = p_st.params;
end




% first component:
x_comp{1} = cos(2 * pi * f0 .* t);
re = zeros(L_freq, L_amp);

% if(strcmp(upper(method), 'EMD'))
%     up_factor = 1;
    
%     % usample for the EMD:
%     x_comp{1} = upsample(x_comp{1}, up_factor);
%     Fs = Fs * up_factor;
% end

if(time_mask)
    mask1 = zeros(1, N);
    mask1(1:floor(N / 3)) = tukeywin(floor(N/3), 0.2);
    x_comp{1} = x_comp{1} .* mask1;
end



for k = 1:L_freq
    for m = 1:L_amp
        fprintf('|k=%d |m=%d | f_ratio=%g | k0=%g | amp_ratio=%g |\n', ...
               k, m, f_ratio(k), k0(k), amp_ratio(m));
        % dispVars(k, m, f_ratio(k), k0(k), amp_ratio(m));
        
        x_comp{2} = amp_ratio(m) .* cos(2 * pi * f1(k) .* t);
        if(time_mask)
            x_comp{2} = x_comp{2} .* fliplr(mask1);
        end


        x = x_comp{1} + x_comp{2};


        if(db_plot)
            err_plot(x, x_comp, x, 9000, 'SIGNALS');            
        end


        try 
            re(k, m) = extract_components(x, Fs, x_comp, method, params, db_plot);
        catch
            re(k, m) = Inf;
        end
        if(db_plot)
          disp('--- paused; hit key to continue ---'); pause;  
        end

    end
end

% set_figure(2000); 
% mesh(f_ratio, amp_ratio, re);

set_figure(2001);
rex = re;
rex(rex > 1) = 1;
imagesc(log10(amp_ratio), f_ratio, rex);
axis('xy');
axis('tight');
xlabel('amplitude ratio (log10)');
ylabel('frequency ratio');

set_figure(2002);
rex = re;
rex(rex > 0.5) = 0.5;
imagesc(log10(amp_ratio), f_ratio, rex);
axis('xy');
axis('tight');
xlabel('amplitude ratio (log10)');
ylabel('frequency ratio');



if(save_)
    if(time_mask)
        tstr = '_septime_tukey';
    else
        tstr = [];
    end

    time_now = now;
    save(['./data/plots/tones_test_' method tstr '_v2.mat'], 're', 'amp_ratio', 'f_ratio', ...
         'time_now');
end





function re = extract_components(x, Fs, x_comp, method, params, db_plot)
%---------------------------------------------------------------------
% decomposition methods and calculate error
%---------------------------------------------------------------------
N_components = length(x_comp);
N = length(x);


[y, y_comps] = decomp_all_methods(x, Fs, method, N_components, params, db_plot);


y = y(:)';
[re, snr] = cal_relative_err(y_comps, x_comp, N_components);    
if(db_plot)
    err_plot(x, y_comps, y, 1000, method);
end






function [re, snr] = cal_relative_err(y_comps, x_components, N_components)
%---------------------------------------------------------------------
% relative error
%---------------------------------------------------------------------
rerr_tmp = [];
snr_tmp = [];

for p = 1:1
    xc = x_components{p};

    r = zeros(1, N_components);
    for q = 1:N_components
        r(q) = corr(y_comps{q}(:), xc(:));
    end
    [~, isort] = max(r);

    diff_l2 = sqrt(sum(abs(xc(:).' - y_comps{isort(1)}(:).') .^ 2));
    comp_l2 = sqrt(sum(abs(x_components{2}(:).') .^ 2));
    rerr_tmp = [rerr_tmp (diff_l2 / comp_l2)];        
    
end
re = mean(rerr_tmp);
snr = 10 * log10(mean(1 ./ rerr_tmp));

print_table([re; snr]', {'RE', 'SNR'});



