function tfrep = estimate_tfd(signal1, wl)
% This function estimates a TFD from a signal - dodgily. There is nothing
% fancy here just a simple separable kernel applied in teh ambiguity
% domain: the window in tau (lag) is twice as wide as in nu (frequency lag)
% to emphasize oscillating signal components. Will replace with Johnno's
% stuff when the time comes. Matlab also has a WVD now as well.
%
% Inputs: signal1 - the time series to be TFDed
%         wl - the window size in the ambguity domain in 2D = wl x 2*wl in
%         span a 2D Hanning window is used
%         
% Outputs: tfrep - the TFD
%
% - dependencies ambig.m
%
% Nathan Stevenson
% July 2021

N = length(signal1);
ref1 = floor(N/2+1);
am2 = hanning(wl)*hanning(2*wl+1)';
am_sep = zeros(N, N);
am_sep(ref1-floor(wl/2):ref1+floor(wl/2), ref1-floor(wl):ref1+floor(wl)) = am2;
amb = ambig(signal1); 
tf_tau = ifft(ifftshift((amb.*am_sep'),2),[],2); 
tfrep = real(1./N.*fft(ifftshift(tf_tau,1), N, 1));