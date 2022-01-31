%-------------------------------------------------------------------------------
% preprocessing_signal: preprocess signal
%
% Syntax: y = preprocessing_signal(x, pad_signal, lpf_signal)
%
% Inputs: 
%     x, pad_signal, lpf_signal - 
%
% Outputs: 
%     y - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 15-09-2021
%
% last update: Time-stamp: <2021-09-15 17:59:38 (otoolej)>
%-------------------------------------------------------------------------------
function y = preprocessing_signal(x, pad_signal, lpf_signal)
if(nargin < 2 || isempty(pad_signal)), pad_signal = false; end
if(nargin < 3 || isempty(lpf_signal)), lpf_signal = false; end

y = x;

if(pad_signal)
    N = length(x);
    L = floor(N / 16);
    
    % positive mirror extrapolation    
    wpad = hamming(2 * L - 1)';
    y = [wpad(1:L) .* fliplr(x(2:L + 1)) x fliplr(wpad(1:L) .* x(end - L:end - 1))];     
    % x = [zeros(1, L) x zeros(1, L)];         
    % x = [x(end - L:end - 1) x x(2:L)];             
end


if(lpf_signal)
    
end

