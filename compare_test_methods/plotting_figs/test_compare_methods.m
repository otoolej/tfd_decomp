% all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'vmd', 'vncmd'};
% [x, x_components, y, y_comps] = compare_methods_testsignals('bat', all_methods, false);
% [x, x_components, y, y_comps] = compare_methods_testsignals('nnlfm4', all_methods, false);
% [x, x_components, y, y_comps] = compare_methods_testsignals('noise', all_methods, false);


% for the plot:
N = 400;
Nq = floor(N / 4);
x1 = gen_LFM(N, 0.05, 0.30);
x2 = gen_LFM(N, 0.25, 0.40, pi / 4);
x1((N - Nq):end) = 0;

x2(1:Nq) = 0;
x = x1 + x2;

[~, ~, Fs, all_params] = set_signal_parameters('bat', false);
params = all_params(2).params;

params.doppler_kernel = {[39], 'hamm'};
params.lag_kernel = {[109], 'dolph', [100]};
[y, y_comps, comp_time] = decomp_all_methods(x, Fs, 'tvfilt', 4, params, true);



