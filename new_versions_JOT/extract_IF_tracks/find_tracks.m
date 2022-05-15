%-------------------------------------------------------------------------------
% find_tracks: wrapper for McAulay--Quatieri method [1] to estimate IF components.
%
% Syntax: [if_tracks,f_scale,t_scale,if_tracks] = find_tracks(tfd,Fs,N,params,freq_limits)
%
% Inputs: 
%     tfd         - TFD matrix
%     Fs          - sampling frequency
%     N           - signal length
%     params      - parameters (class)
%     freq_limits - upper and lower bound frequencies to limit search area 
%
% Outputs: 
%     if_law    - cell of IF tracks 
%     f_scale   - frequency scaling factor
%     t_scale   - time scaling factor
%     if_tracks - unprocessed IF tracks
%     if_energy - energy (sum of TFD) along IF ridge
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
% last update: Time-stamp: <2022-04-20 16:53:49 (otoolej)>
%-------------------------------------------------------------------------------
function [if_law, f_scale, t_scale, if_tracks, if_energy] = ...
    find_tracks(tfd, Fs, N, params, freq_limits)
if(nargin<4 || isempty(params)) params = decomp_params; end
if(nargin<5 || isempty(freq_limits)) freq_limits=[]; end


if_law = []; if_tracks = [];

%---------------------------------------------------------------------
% 0. set parameters
%---------------------------------------------------------------------
[Ntime, Nfreq] = size(tfd);

dur = N / Fs;
t_scale = (dur / Ntime);
f_scale = (1 / Nfreq) * (Fs / 2);

% delta_freq_samples = floor( (params.delta_search_freq / f_scale) * t_scale );
% min_component_length = floor( params.min_if_length / t_scale );
delta_freq_samples = floor(params.delta_search_freq * (Nfreq / Ntime));

% dispVars(params.delta_search_freq, delta_freq_samples, params.min_if_length, Nfreq);


% if are limiting search to specific frequency band:
if(~isempty(freq_limits))
  klimits = ceil( freq_limits(1) / f_scale ):floor( freq_limits(2) / f_scale );

  tfd_narrow = zeros(size(tfd));
  tfd_narrow(:, klimits) = tfd(:, klimits);
  tfd_use = tfd_narrow;
else
  tfd_use = tfd;
end


    
%---------------------------------------------------------------------
% 1. Extract tracks: McAuley--Quatieri method (circa 1986)
%---------------------------------------------------------------------
if_tracks = tracks_IF_MCQ_method(tfd_use, delta_freq_samples, ...
                                 params.min_if_length, params.max_no_peaks);

if(isempty(if_tracks))
    return;
end
    


%---------------------------------------------------------------------
% 2. Order the tracks, in terms of energy (highest energy first)
%    Tracks need to be a certain length.
%---------------------------------------------------------------------
Ntracks = length(if_tracks);

y = []; if_energy = [];

p = 1;
for i = 1:Ntracks
    if(size(if_tracks{i}, 1) > 0)
        ifest = if_tracks{i}(:, 2);
        if_length = length(if_tracks{i}(:, 1));        

        % if between the range we are interested in
        if( if_length >= params.min_if_length )
            n = if_tracks{i}(:, 1);
            k = if_tracks{i}(:, 2);

            subset_tracks{p}(:, 1) = n;
            subset_tracks{p}(:, 2) = k;
            
            % features for each IF:
            if_energy(p) = 0;
            for q = 1:length(n)
                if_energy(p) = if_energy(p) + sum( tfd(n(q), k(q)) );
            end
            
            p = p + 1;
        end
    end
end


% rank IF in relation to energy:
if(~isempty(if_energy))
    [d, imax_energy] = sort(if_energy, 'descend');
    for a = 1:length(imax_energy)
        imaxi = imax_energy(a);

        if(d(a)<0)
            if_law{a} = []; 
        else
            if_law{a} = zeros(length(subset_tracks{imaxi}(:, 1)), 2);    
            if_law{a}(:, 1) = subset_tracks{imaxi}(:, 1);
            if_law{a}(:, 2) = subset_tracks{imaxi}(:, 2);    
        end
    end

    % better to remove IF laws with negative energy:
    if_energy = d;
    irem = find(if_energy <= 0);
    if(~isempty(irem))
        if_energy(irem) = [];
        if_law(irem) = [];
    end

end


