%-------------------------------------------------------------------------------
% extract_components: Extract components from signal using STFT or Q-TFD
%
% Syntax: [y_components]=extract_components(x,Fs,N_components)
%
% Inputs: 
%     x,Fs - 
%
% Outputs: 
%     [y_components] - 
%
% Example:
%     
%

% John M. O' Toole, University of Deusto
% Started: 17-01-2013
%-------------------------------------------------------------------------------
function [y, x_residual, x_components] = ...
    extract_components_xTFD(x, Fs, params, N_components, db_plot, tfd_limits_freq)
if(nargin < 2 || isempty(Fs)) Fs = 1; end
if(nargin < 3 || isempty(params))
    params = decomp_params(length(x), 'xtfd');
end
if(nargin < 4 || isempty(N_components)) N_components = []; end
if(nargin < 5 || isempty(db_plot)) db_plot = 0; end
if(nargin < 6 || isempty(tfd_limits_freq)) tfd_limits_freq = []; end


x = x(:).';

y = []; x_residual = []; x_components = [];
N = length(x);

if(params.pad_signal)
    No = N;
    xo = x;
    x = preprocessing_signal(x, true);
    Lpad = length(x) - N;
    N = length(x);
end



%---------------------------------------------------------------------
% 1. generate the Q-TFD
%---------------------------------------------------------------------
% generate the TFD
[qtfd, g2] = qtfd_sep_kern(x, params.doppler_kernel, params.lag_kernel, N, params.Nfreq);
% scale
qtfd = scale_tfd(qtfd, g2, N);


if(db_plot)
    set_figure(2);
    title('time-frequency distribution');
    vtfd(qtfd); axis('tight');
end


%---------------------------------------------------------------------
% 2. estimate tracks (from peaks) and pick components in TF domain
%---------------------------------------------------------------------
[if_tracks, ~, ~] = find_tracks(qtfd, Fs, N, params, tfd_limits_freq);
% fprintf('n-decomp: %d\n', length(if_tracks));


% limit number of components:
if_tracks = if_tracks(1:min(N_components, length(if_tracks)));


if(isempty(if_tracks))
    if(params.db_warn)
        fprintf('no tracks found in the TFD; exiting\n');
    end

    return;    
end


%---------------------------------------------------------------------
% 3. synthesize components
%---------------------------------------------------------------------
% this turns off correction of the IA from estimates of the BW of the components
if(~params.correct_amplitude_bw)
    g2 = [];
end

% a. estimate components from IF (but without IP):
s_nophase = synth_signal_sinmodel(if_tracks, [], qtfd, Fs, [], []);
%  ... and combine into 1 signal:
% s_all = combine_one_signal(s_nophase);
s_all = sum(cat(2, s_nophase{:}), 2);


% b. X-TFD between original signal and sin-model (without IP):
if(params.phase_correction)

    xtf = xtfd_sep_kern(x, s_all, params.doppler_kernel, params.lag_kernel, N, params.Nfreq);
    
    % c. estimate phase from X-TFD (using IF locations):
    ip_tracks = get_phase(if_tracks, xtf);

    % d. synthesize components again but this time with IPs:
    s_est = synth_signal_sinmodel(if_tracks, ip_tracks, qtfd, Fs, [], g2);

    
    db_testing_plot = false;
    if(db_testing_plot)
        set_figure(987); 
        plot(ip_tracks{1}(:,1), ip_tracks{1}(:,2));
    end

else
    % only for comparing without phase adjustment:
    s_est = s_nophase;
end


% if s_est is not quite the right length:
L_comp = length(s_est);
x_components = cell(L_comp, 1);
for a = 1:length(s_est)
  x_components{a} = zeros(1, N);
  nn = 1:min(N, length(s_est{a}));
  
  x_components{a}(nn) = s_est{a}(nn);
end


%---------------------------------------------------------------------
% 4. remove components from signal
%---------------------------------------------------------------------
if(~isempty(N_components))

    y = nansum(cat(1, x_components{1:min(N_components, length(x_components))}), 1);
    y = y(:)';

    if(length(y) > N)
        y = y(1:N);
    end

    x_residual = x - y;
end



if(params.pad_signal)
    y = y(floor(Lpad / 2) + 1:end - floor(Lpad / 2));
    x_residual = x_residual(floor(Lpad / 2) + 1:end - floor(Lpad / 2));
    for p = 1:length(x_components)
        x_components{p} = x_components{p}(floor(Lpad / 2) + 1:end - floor(Lpad / 2));
    end

    x = xo;
end


%---------------------------------------------------------------------
% 5. PLOT
%---------------------------------------------------------------------
if(db_plot)
    if(params.phase_correction)
        set_figure(301);
        title('Phase of the cross time-frequency distribution');
        vtfd(angle(xtf)); axis('tight');
    end
end



function y = combine_one_signal(x_comps)
%---------------------------------------------------------------------
% combine components into 1 signal
%---------------------------------------------------------------------
y = zeros(length(x_comps{1}), 1);
for a = 1:length(x_comps)
    y = y + x_comps{a};
end

