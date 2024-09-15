function decomp = tf_decomposition_tvfilt(signal1, params, N_components, db_plot)
%
%
%  Note TFD will be NxN (initially) where N is the signal length; TFD,
%  edge linking, and tv-filtering parameters will be in samples and
%  therefore with respect to N
%
%  TFD inputs:
%  wl - window length of separable kernel - 2D size is (wl, 2*wl) to emphasize continuity/smoothness in frequency
%  flag 1 - for some stupid reason I permit an option to whiten the signal
%  before hand, whitening options at this stage is only differentiation.
%  Can change. 0 = don't whiten, 1 = whiten
%
%  FIND COMPONENT inputs:
%  delta_freq_samples - search range in terms of frequency/bandwidth in samples
%  fw - search range in terms of time
%  min_if_length - minimum component length
%
% 
%
if(nargin < 2 || isempty(params))
    params = decomp_params(length(signal1), 'tvfilt');
end
if(nargin < 3 || isempty(N_components)), N_components = []; end
if(nargin < 4 || isempty(db_plot)), db_plot = false; end


N = length(signal1);



%---------------------------------------------------------------------
% zero pad at the end:
%---------------------------------------------------------------------
signal1 = [signal1 zeros(1, params.L_filt + 1)];
N = length(signal1);

%---------------------------------------------------------------------
% generate separable-kernel TFD:
%---------------------------------------------------------------------
[qtfd, g2] = qtfd_sep_kern(signal1, params.doppler_kernel, params.lag_kernel, N, N);
qtfd = scale_tfd(qtfd, g2, N);

db_plot = false;
if(db_plot)
    set_figure(2001); 
    vtfd(qtfd);
    axis('tight');
end

%---------------------------------------------------------------------
% extract IF components from TFD
%---------------------------------------------------------------------
[el1, ei1] = if_tracks_MCQmethod(qtfd, N, 1, ...
                                 params.delta_freq_samples, ...
                                 params.min_if_length, ...
                                 params.max_no_peaks, ...
                                 params.qtfd_max_thres);


%---------------------------------------------------------------------
% do the TV filtering
%---------------------------------------------------------------------
decomp = tv_filtering(el1, ei1, signal1, qtfd', params.L_filt, 'lossless', N_components);


% remove zero-padding
Mh = ceil(params.L_filt / 2);
decomp = decomp(:, Mh:end - Mh - 1);
