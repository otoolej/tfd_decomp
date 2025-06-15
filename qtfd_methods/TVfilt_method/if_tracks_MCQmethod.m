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
% last update: Time-stamp: <2025-06-15 21:56:11 (otoolej)>
%-------------------------------------------------------------------------------
function [if_tracks, tfd_mask] = if_tracks_MCQmethod(qtfd, N, L_filt, Fs, bw, min_if_length, ...
                                                     max_peaks, qtfd_max_thres, freq_limits)
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
if(length(if_tracks) == 0)
    tfd_mask = zeros(size(qtfd, 1), N + L_filt);
    return; 
end



% swap time and frequency columns:
if_tracks = cellfun(@(x) fliplr(x), if_tracks, 'un', false);


%---------------------------------------------------------------------
% 3. ensure that the IF tracks extend cover the whole time course 
%---------------------------------------------------------------------
N_components = length(if_tracks);
[n_min, i_min] = min(cellfun(@(x) min(x(:, 2)), if_tracks));
if n_min > 1;
    cur_if = if_tracks{i_min(1)};
    if_tracks{i_min} = [cur_if(1, 1) * ones(n_min - 1,1) [1:(n_min - 1)]'; cur_if];
end
% ensure that a) IFs extend all the way to then end and b) padded for later filtering
for ii = 1:N_components
    cur_if = if_tracks{ii};
    N_end = max(cur_if(:, 2));
    N_extra = N_end + L_filt;
    if(N_extra > N)
        % if the IF-track needs to padded to allow for the filter convolution:
        if_tracks{ii} = [cur_if; cur_if(end, 1)*ones(L_filt, 1) [N_end+1:N_extra]'];
    end
end
[n_max, i_max] = max(cellfun(@(x) max(x(:, 2)), if_tracks));
N_full = N + L_filt;
if(n_max < N_full)
    cur_if = if_tracks{i_max(1)};
    if_tracks{i_max} = [cur_if; cur_if(end, 1)*ones(N_full - n_max, 1) [n_max+1:N_full]'];
end
%  last, check to see if IFs extend all time:
% i_time = cellfun(@(x) unique(x(:, 2)), if_tracks, 'un', false);
% n_times = sort(unique(vertcat(i_time{:})));
% all_n = 1:N;
% n_missing = all_n(~ismember(all_n, n_times));
% if(length(n_missing))
%     new_if = [floor(N / 2) * ones(length(n_missing)) n_missing];
%     if_tracks{length(if_tracks) + 1} = new_if;
% end



%---------------------------------------------------------------------
% 4. generate a TFD mask along the IF ridges
%---------------------------------------------------------------------
freq_maxes = cellfun(@(x) max(x(:, 1)), if_tracks);
tfd_mask = zeros(size(qtfd, 1), N_full);
for n = 1:N_components
    cur_if = if_tracks{n};
    tfd_mask(sub2ind(size(tfd_mask), cur_if(:, 1), cur_if(:, 2))) = 1; 
end


dbplot = false;
if(dbplot)
    set_figure(2002); 
    mesh(tfd_mask);
end


