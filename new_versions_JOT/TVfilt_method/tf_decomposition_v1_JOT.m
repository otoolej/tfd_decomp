function decomp = tf_decomposition_v1_JOT(signal1, params, N_components, db_plot)
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
%  bw - search range in terms of frequency/bandwidth in samples
%  fw - search range in terms of time
%  len1 - minimum component length
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
qtfd = qtfd_sep_kern(signal1, params.doppler_kernel, params.lag_kernel, N, N);
% qtfd = qtfd_sep_kern(signal1, {params.wx, 'hann'}, {params.wx, 'dolph', 100}, N, N);


% dispVars(size(tfrep));
if(db_plot)
    set_figure(2001); 
    vtfd(qtfd);
    axis('tight');
end


%---------------------------------------------------------------------
% extract IF components from TFD
%---------------------------------------------------------------------
% [el1, ei1] = find_components_v2(qtfd', qtfd', params.bw, params.fw, params.len1);
[el1, ei1] = if_tracks_MCQmethod(qtfd, N, 1, ...
                                 params.bw, params.fw, params.len1, ...
                                 params.max_no_peaks, ...
                                 params.qtfd_max_thres);


%---------------------------------------------------------------------
% do the TV filtering
%---------------------------------------------------------------------
decomp = tv_filtering_JOT(el1, ei1, signal1, qtfd', params.L_filt, 'lossless', N_components);


% remove zero-padding
Mh = ceil(params.L_filt / 2);
decomp = decomp(:, Mh:end - Mh - 1);
