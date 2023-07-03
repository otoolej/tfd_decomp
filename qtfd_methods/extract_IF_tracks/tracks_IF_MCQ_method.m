%-------------------------------------------------------------------------------
% tracks_IF_MCQ_method: wrapper for McAuley-Quatieri method
%
% Syntax: tf_tracks = tracks_IF_MCQ_method(tf, delta_limit, min_length, max_npeaks)
%
% Inputs: 
%     tf          - TFD matrix
%     delta_limit - max. rate of change of IF track 
%     min_length  - minimum length of IF track
%     max_npeaks  - max. number of peaks to estimate per time-slice of the TFD
%
% Outputs: 
%     individual_tracks - IF tracks (cell)
%
% Example:
%     
%
% [1] RJ McAulay, TF Quatieri. Speech analysis/synthesis based on a sinusoidal
% representation. IEEE Transactions on Acoustics, Speech, and Signal Processing, 34(4),
% (1986), 744–754.


% John M. O' Toole, University College Cork
% Started: 27-08-2021
%
% last update: Time-stamp: <2021-08-27 13:05:27 (otoolej)>
%-------------------------------------------------------------------------------
function individual_tracks = tracks_IF_MCQ_method(tf, delta_limit, min_length, max_npeaks)
if(nargin < 2 || isempty(delta_limit)) delta_limit = 4; end
if(nargin < 3 || isempty(min_length)) min_length = []; end
if(nargin < 4 || isempty(max_npeaks)), max_npeaks = 8; end


%---------------------------------------------------------------------
% estimate the IFs and return as 'tracks'
%---------------------------------------------------------------------
mq_params.MAX_NO_PEAKS = max_npeaks; % max. number of peaks to extract
mq_params.delta = delta_limit;  % and max. rate of change of IF 
[~, individual_tracks] = extract_IF_tracks(tf, mq_params);

%---------------------------------------------------------------------
% 2.A) MODIFY: remove tracks of length less than three
%      as with birth and death this is only really a 1-point peak
%---------------------------------------------------------------------
if(~isempty(min_length) & isempty(individual_tracks))
    N_tracks = length(individual_tracks);
    nmat = [];
    for n = 1:N_tracks
        if(length(individual_tracks{n}) <= min_length)
            nmat = [nmat, n];
        end
    end
    individual_tracks(nmat) = [];
end
