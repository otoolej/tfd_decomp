function decomp = tf_decomposition_v1(signal1, len1)
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

N = length(signal1);
[wx, bw, fw, M] = setparameters(N);
bwf = 'lossless';

% estimate TFD
%switch flag1
%    case 0
      tfrep = estimate_tfd(signal1, wx);
      tfrep1 = tfrep;
%    case 1
%       tfrep = estimate_tfd(signal1(1:N-1), wx);
%       tfrep1 = estimate_tfd(diff(signal1), wx); 
%end

% find components in TFD
[el1, ei1] = find_components(tfrep, tfrep1, bw, fw, len1);

% extract components from signal
decomp = tv_filtering(el1, ei1, signal1, tfrep, M, bwf);
decomp = decomp(:, ceil(M/2):end);
