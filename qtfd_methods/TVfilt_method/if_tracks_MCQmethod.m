%-------------------------------------------------------------------------------
% if_tracks_MCQmethod: wrapper to find_tracks.m file for TV-filt method
%
% Syntax: [if_tracks, tfd_mask] = if_tracks_MCQmethod(qtfd, N, Fs, bw, min_if_length, max_peaks, freq_limits)
%
% Inputs: 
%     qtfd        - TFD matrix
%     N           - signal length
%     Fs          - sampling frequency
%     bw          - maximum frequency shift to allow when joining tracks (in samples)
%     min_if_length  - minimum length of IF track (in samples)
%     max_peaks   - maximum number of peaks per time slice
%     freq_limits - upper and lower bound frequencies to limit search area 
% 
%
% Outputs:
%     if_tracks   - cell of IF tracks 
%     tfd_mask    - binary matrix with TF location of IF tracks
%
% Example:
%     
%
% [1] RJ McAulay, TF Quatieri. Speech analysis/synthesis based on a sinusoidal
% representation. IEEE Transactions on Acoustics, Speech, and Signal Processing, 34(4),
% (1986), 744–754.


% John M. O' Toole, University College Cork
% Started: 14-04-2022
%
% last update: Time-stamp: <2023-07-19 19:25:45 (otoolej)>
%-------------------------------------------------------------------------------
function [if_tracks, tfd_mask] = if_tracks_MCQmethod(qtfd, N, Fs, bw, min_if_length, max_peaks, ...
                                                  qtfd_max_thres, freq_limits)
if(nargin < 3 || isempty(Fs)), Fs = 1; end
if(nargin < 8 || isempty(qtfd_max_thres)), qtfd_max_thres = []; end
if(nargin < 9 || isempty(freq_limits)), freq_limits = []; end

%---------------------------------------------------------------------
% 0. remove low-amplitude activity in the TFD
%---------------------------------------------------------------------
if(~isempty(qtfd_max_thres))
    th = qtfd_max_thres * max(max(qtfd));
    % a = size(qtfd);
    % th = sum(sum(qtfd))/(2*prod(a));
    % th = mean(qtfd(:)) / 2;
    qtfd(qtfd < th) = 0;
end



%---------------------------------------------------------------------
% 1. set the parameters
%---------------------------------------------------------------------
params_mcq.delta_freq_samples = bw;
params_mcq.min_if_length = min_if_length;
params_mcq.max_no_peaks = max_peaks;

%---------------------------------------------------------------------
% 2. extract the IF components
%---------------------------------------------------------------------
[if_tracks, ~, ~, ~, if_energy] = find_tracks(qtfd, Fs, N, params_mcq, freq_limits);

% swap time and frequency columns:
if_tracks = cellfun(@(x) fliplr(x), if_tracks, 'un', false);


%---------------------------------------------------------------------
% 3. ensure that the IF tracks extend cover the whole time course 
%---------------------------------------------------------------------
m1 = zeros(1,length(if_tracks));
m2 = zeros(1,length(if_tracks));
for ii = 1:length(if_tracks)
   m1(ii) = min(if_tracks{ii}(:, 2));
   m2(ii) = max(if_tracks{ii}(:, 2));
end
if m1>1;
    aa = find(m1==min(m1));
    aa = aa(1);
    if_tracks{aa} = [if_tracks{aa}(1,1)*ones(min(m1),1) [1:min(m1)]' ; if_tracks{aa}];
end
if m2<N;
    aa = find(m2==max(m2));
    aa = aa(1);
    if_tracks{aa} = [if_tracks{aa} ; if_tracks{aa}(end,1)*ones(N-max(m2),1) [max(m2)+1:N]'];
end


%---------------------------------------------------------------------
% 4. generate a TFD mask along the IF ridges
%---------------------------------------------------------------------
tfd_mask = zeros(size(qtfd));

for n = 1:length(if_tracks)
    if_ = if_tracks{n};
    for p = 1:length(if_)
        tfd_mask(if_(p, 1), if_(p, 2)) = 1;
    end
end


dbplot = false;
if(dbplot)
    set_figure(2002); 
    mesh(tfd_mask);
end

    
