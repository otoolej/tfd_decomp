% STENNY BASED DEMO
% JUST WORKING ON TWO EXAMPLES SIGNALS WHICH I AM PUTTING AS EXTREMES
% A COUPLE OF COMPONENT - BAT and A LOT OF COMPONENTS - NOISE
%cd('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JohnO\tfd_decomp-master')
%addpath(genpath('L:\Lab_JamesR\nathanST\QIMRBerghofer_Stevenson\Document_and_Code\people\Australia\QIMRB\JohnO\tfd_decomp-master'))

d = load('batsignal');
x = d.batsignal; 
N = length(d.batsignal)
params = decomp_params(length(x), 'XTFD');
params.max_no_peaks = Inf;
params.lag_kernel = {63, 'dolph', 100};
params.doppler_kernel = {63 'hamm'};    
[qtfd, g2] = qtfd_sep_kern(x, params.doppler_kernel, params.lag_kernel, N, params.Nfreq);
[qtfd, scale_factor] = scale_tfd(qtfd, g2, N);
[if_tracks, ~, ~] = find_tracks(qtfd, 1./d.DT, N, params, []);
% problematic methinks - way too many tracks selected so we need some sort of
% prune - easiest way is to hack the qtfd for starters so dump negative
% energy or even better dump at a threshold related to auniform spread of energy
% i.e set a threshold equivalent to 50% of the sum(sum(qtfd))./prod(size(qtfd))
% this should also speed up track estimation
figure; hold on;
for ii = 1:length(if_tracks)
    plot(if_tracks{ii}(:,1), if_tracks{ii}(:,2))
end
th = sum(sum(qtfd))./prod(size(qtfd))/2;
length(if_tracks)
qtfd(qtfd<th)=0;
[if_tracks, f_scale, t_scale] = find_tracks(qtfd, 1./d.DT, N, params, []);
figure; hold on;
for ii = 1:length(if_tracks)
    plot(if_tracks{ii}(:,1), if_tracks{ii}(:,2))
end
% There should also be some sort of limit in duration of tracks, purely
% related to smoothing window size, i.e. if I have a TF delta and do the 2D
% smoothing then there is a limit lower limit in bandwidth and duration
bw = zeros(1,length(if_tracks)); dur = zeros(1,length(if_tracks));
for ii = 1:length(if_tracks)
    bw(ii) = max(if_tracks{ii}(:,1)) - min(if_tracks{ii}(:,1));
    dur(ii) = max(if_tracks{ii}(:,2)) - min(if_tracks{ii}(:,2));
end
th = sqrt(size(qtfd));
ref = find(bw>th(1) & dur>th(2));
figure; hold on;
for ii = 1:length(ref)
    plot(if_tracks{ref(ii)}(:,1), if_tracks{ref(ii)}(:,2))
end
% considering how much track you are interested in you can then go and find
% where these tracks overlap the original tracks and use these, I am not
% sure if the xTFD need longer or shorter tracks.
% NOTE THAT ALL THIS FAFFIN ABOUT SHOULD LIMIT THE NUMBER OF TRACKS TO an
% upper limit of prod(th) - if I assume an N x N TFD, the limit will be
% N which makes sense

% Let's see if the same procedure works for the fractional Gaussian noise
% case as this on and the bat signal would seem to be the extreme
load('noisysignals.mat');
params = decomp_params(length(x), 'XTFD');
params.pad_signal = 0;
params.max_no_peaks = Inf;
params.lag_kernel = {63, 'dolph', 100};
params.doppler_kernel = {63, 'hamm'};
params.delta_search_freq = 1000;
N = length(x);

[qtfd, g2] = qtfd_sep_kern(x, params.doppler_kernel, params.lag_kernel, N, params.Nfreq);
[qtfd, scale_factor] = scale_tfd(qtfd, g2, N);
[if_tracks, ~, ~] = find_tracks(qtfd, 1, N, params, []);
% problematic methinks - way too many tracks selected so we need some sort of
% prune - easiest way is to hack the qtfd for starters so dump negative
% energy or even better dump at a threshold related to uniform spread of energy
% i.e set a threshold equivalent to 50% of the sum(sum(qtfd))./prod(size(qtfd))
% this should also speed up track estimation
figure; hold on;
for ii = 1:length(if_tracks)
    plot(if_tracks{ii}(:,1), if_tracks{ii}(:,2))
end
th = sum(sum(qtfd))./prod(size(qtfd))/2; % new threshold
length(if_tracks)
qtfd(qtfd<th)=0;
[if_tracks, f_scale, t_scale] = find_tracks(qtfd, 1, N, params, []);
figure; hold on;
for ii = 1:length(if_tracks)
    plot(if_tracks{ii}(:,1), if_tracks{ii}(:,2))
end
% Also need to do some scanning to prevent disconinuity in IF laws of which
% seem to happen in the random case so added this bit, once again I will be
% guided by sqrt(N) here - seems to work and not effect deterministic signals 
th = ceil(sqrt(size(qtfd))); M = length(if_tracks);
for ii = 1:length(if_tracks)
    cc = if_tracks{ii}(:,2);
    dd = diff(cc);
    if max(abs(dd))>2*th(2) % then split track
        rx = find(abs(dd)==max(abs(dd)));
        if_tracks{ii} = if_tracks{ii}(1:rx,:);
        if_tracks{M+1} = if_tracks{ii}(rx+1:end,:);
        M = M+1;
     end
end
% There should also be some sort of limit in duration of tracks, purely
% related to smoothing window size, i.e. if I have a TF delta and do the 2D
% smoothing then there is a limit lower limit in bandwidth and duration
% in this case, it is too much hard work so I have just assume eqal
% smoothing in both t and f directions
bw = zeros(1,length(if_tracks)); dur = zeros(1,length(if_tracks));
for ii = 1:length(if_tracks)
    bw(ii) = max(if_tracks{ii}(:,1)) - min(if_tracks{ii}(:,1));
    dur(ii) = max(if_tracks{ii}(:,2)) - min(if_tracks{ii}(:,2));
end
th = sqrt(size(qtfd));
ref = find(bw>th(1) & dur>th(2));
figure; hold on;
for ii = 1:length(ref)
    plot(if_tracks{ref(ii)}(:,1), if_tracks{ref(ii)}(:,2))
end



