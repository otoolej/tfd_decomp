%-------------------------------------------------------------------------------
% synth_sumofsins: Synthesis signal from IF information (tracks) and TFD (S)
%
% Different to 'synthesize_signal_if_interp.m' as here no interpolation between
% discrete time slices
%
%
% Syntax: s_est=synth_sumofsins(tracks,S)
%
% Inputs: 
%     tracks,S - 
%
% Outputs: 
%     s_est - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 18-09-2020
%
% last update: Time-stamp: <2024-11-02 19:18:43 (otoolej)>
%-------------------------------------------------------------------------------
function s_est_tracks = synth_signal_sinmodel(if_tracks, ip_tracks, S, Fs, ia_tracks, g2)
if(nargin < 2 || isempty(ip_tracks)), ip_tracks=[]; end
if(nargin < 5 || isempty(ia_tracks)), ia_tracks = []; end
if(nargin < 6 || isempty(g2)), g2 = []; end




[Ntime, Nfreq] = size(S);

if(~isreal(S))
    Smag = abs(S);
else
    Smag = S;
end
f_scale = 1 / (2 * Nfreq / Fs);

Ntracks = length(if_tracks);
s_est = zeros(Ntime, 1);


%---------------------------------------------------------------------
% Method assumes that TFD is sampled at correct rate
% (i.e. no interpolation needed)
%---------------------------------------------------------------------
for it = 1:Ntracks
    cur_track = if_tracks{it};
    L = length(cur_track);

    s_est_tracks{it} = zeros(max(Ntime, L), 1);
    % dispVars(it);


    

    %---------------------------------------------------------------------
    % Store the instantaneous amplitude, frequency, and phase:
    %---------------------------------------------------------------------
    amp_track = zeros(L,1); phase_track = zeros(L, 1);    

    % a. instantaneous frequency (IF)
    if_track = cur_track(:, 2) .* f_scale;     

    % b. estimate instantaneous amplitude (IA)
    for mm = 1:L
        amp_track(mm) = Smag(cur_track(mm, 1), cur_track(mm, 2));
    end
    
    % c. if adding phase correction from the X-TFD:
    if(~isempty(ip_tracks))
        phase_track = ip_tracks{it}(:, 2);
    end
    
    %---------------------------------------------------------------------
    % IA adjustment:
    %---------------------------------------------------------------------
    if(~isempty(g2))

        % -3 dB point for of the lag window:
        half_amp = 10 ^ (-3 / 10);
        
        bwq = zeros(L, 1);
        ikeep = [];
        for mm = 1:L
            bwq(mm) = get_bandwidth(Smag(cur_track(mm, 1), :), cur_track(mm, 2), half_amp);
            if(~isempty(bwq(mm)) && bwq(mm) > 1)
                %  bandwidth shouldn't be less than 1 sample point:
                ikeep = [ikeep mm];
            end
        end
        max_amp = max(amp_track);
        bwq_max = bwq(find(amp_track == max_amp, 1));
        bwq = bwq ./ bwq_max; 
        amp_track(ikeep) = amp_track(ikeep) .* bwq(ikeep);
        amp_track(amp_track > max_amp) = max_amp;
        
    end
    amp_track(amp_track < 0) = 0;
    
    
    DBplot_tracks = false;
    if(DBplot_tracks)
        if(it == 1)
            set_figure(444); 
            plot(cur_track(:, 1), amp_track);
            pk_t = amp_track(ceil(L / 2));
            w_bart = bartlett(L);
            
            plot(cur_track(:, 1), ((w_bart) .* pk_t));
            plot(cur_track(:, 1), pk_t .* amp_track .^ (1 / L));
        end
    end

    %---------------------------------------------------------------------
    % assemble the sinusoidal component:
    %---------------------------------------------------------------------
    nn = cur_track(1, 1):cur_track(end, 1);
    % nn = cur_track(:, 1);
    phase_est = 2.0 * pi .* cumtrapz(if_track.') ./ Fs + phase_track.';

    s_est_tracks{it}(nn) = sqrt(amp_track) .* cos(phase_est.');
    
    % in case it runs off the edge...
    s_est_tracks{it} = s_est_tracks{it}(1:Ntime);
    
    s_est = s_est + s_est_tracks{it};
end



