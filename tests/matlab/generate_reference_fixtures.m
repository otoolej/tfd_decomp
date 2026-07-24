% Generate compact golden-reference data for the Python parity tests.
this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(this_file)));
addpath(genpath(repo_root));

fixture_dir = fullfile(repo_root, 'tests', 'fixtures');
if ~exist(fixture_dir, 'dir')
    mkdir(fixture_dir);
end

N = 64;
n = 0:N - 1;
x = cos(2*pi.*(0.06.*n + ((0.24-0.06)/(2*N)).*(n.^2)));

params_xtfd = decomp_params(N, 'xtfd');
params_xtfd.Nfreq = 128;
params_xtfd.min_if_length = 8;

[analytic_z, ~, ~] = get_analytic_signal(x);
[qtfd_raw, g2, G1] = qtfd_sep_kern(...
    x, params_xtfd.doppler_kernel, params_xtfd.lag_kernel, N, params_xtfd.Nfreq);
qtfd_scaled = scale_tfd(qtfd_raw, g2, N);
[tracks, ~, ~] = find_tracks(qtfd_scaled, 1, N, params_xtfd);
[components_no_phase] = synth_signal_sinmodel(tracks, [], qtfd_scaled, 1, [], []);
signal_no_phase = sum(cat(2, components_no_phase{:}), 2);
[cross_distribution] = xtfd_sep_kern(...
    x, signal_no_phase, params_xtfd.doppler_kernel, params_xtfd.lag_kernel, ...
    N, params_xtfd.Nfreq);
[phase_tracks] = get_phase(tracks, cross_distribution);
[components_with_phase] = synth_signal_sinmodel(...
    tracks, phase_tracks, qtfd_scaled, 1, [], g2);
[y_xtfd, components_xtfd] = tfd_decomposition(x, 'xtfd', 3, params_xtfd, false);

params_tvfilt = decomp_params(N, 'tvfilt');
params_tvfilt.min_if_length = 8;
[y_tvfilt, components_tvfilt] = tfd_decomposition(x, 'tvfilt', 3, params_tvfilt, false);

fixture_version = 1;
save(fullfile(fixture_dir, 'matlab_reference.mat'), ...
     'fixture_version', 'x', 'analytic_z', 'g2', 'G1', 'qtfd_raw', ...
     'qtfd_scaled', 'tracks', 'components_no_phase', 'signal_no_phase', ...
     'cross_distribution', 'phase_tracks', 'components_with_phase', ...
     'y_tvfilt', 'components_tvfilt', 'y_xtfd', 'components_xtfd', '-v7');

fprintf('Wrote %s\n', fullfile(fixture_dir, 'matlab_reference.mat'));

% Final-output references using compact segments of every bundled signal type.
data_dir = fullfile(repo_root, 'data', 'test_signals');
white_data = load(fullfile(data_dir, 'white_noise.mat'));
ffgn_data = load(fullfile(data_dir, 'ffgn_1_02_1_512_0_signal.mat'));
eeg_data = load(fullfile(data_dir, 'eeg_signal.mat'));
eeg1_data = load(fullfile(data_dir, 'test_signal_eeg1.mat'));

dataset_signals.white_noise = white_data.x;
dataset_signals.ffgn = ffgn_data.x;
dataset_signals.eeg = eeg_data.eeg{1}(1:512);
dataset_signals.eeg1 = eeg1_data.eeg1(1:512);
bat = load(fullfile(data_dir, 'bat1'));
dataset_signals.bat = bat(:).';
whale = load(fullfile(data_dir, 'whale1'));
dataset_signals.whale = whale(1:512).';

dataset_names = fieldnames(dataset_signals);
dataset_references = struct();
for dataset_index = 1:length(dataset_names)
    dataset_name = dataset_names{dataset_index};
    dataset_x = dataset_signals.(dataset_name);
    dataset_x = dataset_x(1:end - mod(length(dataset_x), 2));

    dataset_tv_params = decomp_params(length(dataset_x), 'tvfilt');
    [dataset_y_tv, dataset_components_tv] = tfd_decomposition(...
        dataset_x, 'tvfilt', 3, dataset_tv_params, false);

    dataset_xtfd_params = decomp_params(length(dataset_x), 'xtfd');
    dataset_xtfd_params.Nfreq = 1024;
    [dataset_y_xtfd, dataset_components_xtfd] = tfd_decomposition(...
        dataset_x, 'xtfd', 3, dataset_xtfd_params, false);

    dataset_references.(dataset_name).x = dataset_x;
    dataset_references.(dataset_name).y_tvfilt = dataset_y_tv;
    dataset_references.(dataset_name).components_tvfilt = dataset_components_tv;
    dataset_references.(dataset_name).y_xtfd = dataset_y_xtfd;
    dataset_references.(dataset_name).components_xtfd = dataset_components_xtfd;
end

save(fullfile(fixture_dir, 'matlab_dataset_reference.mat'), ...
     'fixture_version', 'dataset_references', '-v7');
fprintf('Wrote %s\n', fullfile(fixture_dir, 'matlab_dataset_reference.mat'));
