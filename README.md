# TFD decomposition method

Code for 2 decomposition methods, both of which estimate the instantaneous frequency of
signal components using quadratic time-frequency distributions (TFDs). The first method
decomposes the signal using time-varying filtering (the `TV-filt` method). The second
method decomposes the signal using phase information from the cross TFD (the `xTFD`
method). See reference O'Toole and Stevenson (2025) for full details.

Use the `tfd_decomposition` to call either method. See `help tfd_decomposition` for more
details. For example, to extract 2 components from signal `x` using the TV-filter method, then:

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


% perform the decomposition
[y, y_comps] = tfd_decomposition(x, 'tvfilt', 2);
```

where `y_comps` is a cell array of the 2 components and `y` is sum of the components.


For example, compare the true and estimated first component:
```matlab
figure(1); hold all;
plot(x1);
plot(y_comps{1});				
```


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

---

## Reference 

JM O'Toole and NJ Stevenson, /Nonstationary Signal Decomposition using Quadratic
Time--Frequency Distributions/, Signal Processing, 2025.
