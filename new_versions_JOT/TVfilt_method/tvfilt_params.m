function params = tvfilt_params(N)

M = sqrt(N);

% length of Doppler and lag window:
params.wx = set_odd(M);

% edge-linker parameters (max. bandwidth to search and ?)
params.bw = set_odd(N/params.wx/2);
params.fw = set_odd(N/params.wx/4);
% minimum length of component:
params.len1 = floor(N / 8);
% params.len1 = sqrt(N);

% zero-padding of the signal:
params.M = set_odd(ceil(2*M));



function x = set_odd(x)
x = x + (1 - rem(x, 2));
