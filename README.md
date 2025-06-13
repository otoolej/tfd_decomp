# TFD decomposition method

Code for 2 decomposition methods, both of which estimate the instantaneous frequency of
signal components using quadratic time-frequency distributions (TFDs). The first method
decomposes the signal using time-varying filtering (the `TV-filt` method). The second
method decomposes the signal using phase information from the cross TFD (the `xTFD`
method). See reference O'Toole and Stevenson (2025) for full details.

Use the `tfd_decomposition` to call either method. See `help tfd_decomposition` for more
details. (First add the path for the library, e.g. `addpath(genpath(pwd));`) For example,
to extract 2 components from signal `x` using the TV-filter method, then:

```matlab
%  generate test signal with 2 components
N = 400;
n = 0:N - 1;
Nq = floor(N / 4);
x1 = cos( 2*pi.*(0.05.*n + ((0.3-0.05)/(2*N)).*(n.^2)) );
x2 = cos( 2*pi.*(0.25.*n + ((0.4-0.25)/(2*N)).*(n.^2)) + (pi / 4));
x1((N - Nq):end) = 0;
x2(1:Nq) = 0;
x = x1 + x2;


% perform the decomposition and plot:
[y, y_comps] = tfd_decomposition(x, 'tvfilt', [], [], true);
```

where `y_comps` is a cell array of the 2 components and `y` is sum of the components. The
input arguments for `tfd_decomposition` are as follows:

```matlab
% Syntax: [y, y_comps] = tfd_decomposition(x, method, N_components, params, db_plot)
%
% Inputs: 
%     x            - signal (length-N)
%     method       - decomposition method; either 'xtfd' or 'tvfilt' (default)
%     N_components - maximum number of components to extract (default=10)
%     params       - parameters for the methods (object; see 'decomp_params.m')
%     db_plot      - plot (default=false)
% 
%
% Outputs: 
%     y       - output signal consisting of sum of components
%     y_comps - cell array containing the components
```

Note that the parameter `N_components` controls the maximum number of components
returned. The number of components generated depends on the number of instantaneous
frequency laws extracted from the TFD.


## Parameters
Each method has a set of parameters that can be changed. These are managed through the
class `decomp_params`. To set the parameters, create an `decomp_params` object
passing in both signal length and method, e.g.

```matlab
params_tvfilt = decomp_params(length(x), 'tvfilt');
```
which generates the parameters for the `TV-filt` method for signal `x`. To change the
TFD kernel parameters, for example, do

```matlab
params_tvfilt.doppler_kernel = {151, 'hann'};
params_tvfilt.lag_kernel = {51, 'hamm'};
```
which sets the Doppler kernel to Hanning window of length 151 samples and the lag kernel
to a Hamming window of length 51 samples.

To change the resolution of sampling in the frequency direction in the TFD for the `xTFD` method, do as follows:
	
```matlab
% first generate the parameter set:
params_xtfd = decomp_params(256, 'xtfd'); 

% set the rate of sampling in the frequency direction of the TFD:
params_xtfd.Nfreq = 32768;
```
which sets the TFD to be of size `256 x 32,768`. Accuracy for the xTFD method may improve
by increasing `Nfreq`: a larger value reduces error but increase computational load and memory
requirements. The default value for `Nfreq` is set to 8,192.

### Parameter list

The full list of default parameters for signal of length `N` 

| parameter              | description                                            | default               | specific to method |
|:-----------------------|--------------------------------------------------------|-----------------------|--------------------|
| `max_no_peaks`         | maximum number of peaks to consdier at each time-slice | 8                     |                    |
| `delta_freq_samples`   | maximum rate-of-change of IF                           | √N/2                  |                    |
| `min_if_length`        | minmum length of IF                                    | 4√N                   |                    |
| `doppler_kernel`       | doppler window (TFD kernel)                            | {4√N, 'hamm'}         |                    |
| `lag_kernel`           | lag window (TFD kernel)                                | {2√N, 'dolph', [100]} |                    |
| `pad_signal`           | pad signal before processing                           | false                 |                    |
| `low_pass_filter`      | low-pass filter before processing                      | false                 |                    |
| `correct_amplitude_bw` | apply amplitude correction                             | true                  | XTFD               |
| `phase_correction`     | apply phase correction                                 | true                  | XTFD               |
| Nfreq                  | sampling rate in frequency direction of TFD            | 8192                  | XTFD               |
| `L_filt`               | filter length                                          | 2√N                   | TV-FILT            |
| `qtfd_max_thres`       | threshold TFD plane with this fraction of maximum      | 0.01                  | TV-FILT            |


### Examples of different TFD kernels
The TFD kernel controls the level of smoothing in the TFD and this can directly impact on
the number and type of components extracted. Here are some examples:


```matlab
x = randn(1,1024);
err = zeros(4, 2);
n_max_components = 200;

params_tvfilt = decomp_params(length(x), 'tvfilt');
params_xtfd = decomp_params(length(x), 'xtfd');

% default parameters for both methods produces 8 components:
[r1, r_comps1] = tfd_decomposition(x, 'tvfilt', n_max_components, params_tvfilt);
err(1, 1) = sum(abs(x-r1'))./sum(abs(x));
err(1, 2) = length(r_comps1);

[r1, r_comps1] = tfd_decomposition(x, 'xtfd', n_max_components, params_xtfd);
err(2, 1) = sum(abs(x-r1))./sum(abs(x));
err(2, 2) = length(r_comps1);

% shorter minimum-duration of components generates more components:
params_tvfilt.min_if_length = 32;
[r1, r_comps1] = tfd_decomposition(x, 'tvfilt', n_max_components, params_tvfilt);
err(3, 1) = sum(abs(x-r1'))./sum(abs(x));
err(3, 2) = length(r_comps1);

params_tvfilt.min_if_length = 12;
[r1, r_comps1] = tfd_decomposition(x, 'tvfilt', n_max_components, params_tvfilt);
err(4, 1) = sum(abs(x-r1'))./sum(abs(x));
err(4, 2) = length(r_comps1);

print_table(err, {'error', '#components'}, {'TVFILT', 'XTFD', 'TVFILT (32)', 'TVFILT (12)'}, [], [], [3, 0]);
```



---

## Reference 

JM O'Toole and NJ Stevenson, "Nonstationary Signal Decomposition using Quadratic
Time--Frequency Distributions", Signal Processing, 2025.
