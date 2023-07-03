%-------------------------------------------------------------------------------
% get_phase: Extract phase from X-TFD using IF 'tracks'
%
% Syntax: ip_tracks = get_phase(if_tracks,xtf,Fs)
%
% Inputs: 
%     if_tracks   - IF tracks (cell)
%     xtf         - cross TFD
%     Fs          - sampling frequency
%     mag_instead - estimate magnitude (default = false)
%
% Outputs: 
%     ip_tracks - instantaneous phase tracks (cell)
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 04-10-2015
%-------------------------------------------------------------------------------
function ip_tracks=get_phase(if_tracks, xtf, Fs, mag_instead)
if(nargin < 3 || isempty(Fs)) Fs=1; end
if(nargin < 4 || isempty(mag_instead)), mag_instead = false; end



[Ntime,Nfreq] = size(xtf);
Ntracks = length(if_tracks);

if(mag_instead)
    angle_or_mag = @(x) sqrt(abs(x));
else
    angle_or_mag = @(x) angle(x);
end


for it = 1:Ntracks
    
    cur_track = if_tracks{it};
    L = length(cur_track);
    
    %---------------------------------------------------------------------
    % Store the amp, IF, and phase points:
    %---------------------------------------------------------------------
    cur_ph_track = zeros(L,2);

    cur_ph_track(:,1) = cur_track(1,1):cur_track(end,1);
    % cur_ph_track(:,1) = cur_track(:, 1);
    for mm = 1:L
        if(mm <= Ntime & cur_track(mm,1) <= Ntime)
            cur_ph_track(mm, 2) = angle_or_mag( xtf(cur_track(mm, 1), cur_track(mm, 2)) );
        end
    end

    ip_tracks{it} = cur_ph_track;
end




