function params = tvfilt_params(N)

M = sqrt(N);

% length of Doppler and lag window:
params.wx = make_odd(M);

% edge-linker parameters (max. bandwidth to search and ?)
params.bw = make_odd(N/params.wx/2);
params.fw = make_odd(N/params.wx/4);
% minimum length of component:
params.len1 = floor(N / 8);
% params.len1 = sqrt(N);

% zero-padding of the signal:
params.M = make_odd(ceil(2*M));


