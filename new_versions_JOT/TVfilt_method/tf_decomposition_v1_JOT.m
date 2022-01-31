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
%  EXTRACT COMPONENT inputs:
%  M = filter order in samples
%  flag 2 - lossy or lossless decomposition (0 = lossy, 1 = lossless)
% 
%
if(nargin < 2 || isempty(params)), params = tvfilt_params(length(signal1)); end

N = length(signal1);


% zero pad at the end:
signal1 = [signal1 zeros(1, params.M + 1)];
N = length(signal1);

% generate separable-kernel TFD:
Ldopp = make_odd(floor(1.5 * params.wx));
qtfd = qtfd_sep_kern(signal1, {Ldopp, 'hann'}, {params.wx, 'dolph', 100}, N, N);
% qtfd = qtfd_sep_kern(signal1, {params.wx, 'hann'}, {params.wx, 'dolph', 100}, N, N);
tfrep = qtfd';


% dispVars(size(tfrep));
if(db_plot)
    set_figure(2001); 
    vtfd(qtfd);
    axis('tight');
end

      
% find components in TFD
[el1, ei1] = find_components_v2(tfrep, tfrep, params.bw, params.fw, params.len1);

% extract components from signal
decomp = tv_filtering_JOT(el1, ei1, signal1, tfrep, params.M, 'lossless', N_components);

% remove zero-padding
Mh = ceil(params.M / 2);
decomp = decomp(:, Mh:end - Mh - 1);
