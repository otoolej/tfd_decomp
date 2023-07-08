# TFD decomposition method

Code for 2 decomposition methods, both of which estimate the instantaneous frequency of
signal components using quadratic time-frequency distributions (TFDs). The first method
decomposes the signal using a time-varying filtering (the `TV-filt` method). The second
method decomposes the signal using phase information from the cross TFD (the `xTFD`
method).

Use the `tfd_decomposition` to call either method. See `help tfd_decomposition` for more
details. For example, to extract 2 components from signal `x` using the TV-filter method, then:

```matlab
[y, y_comps] = tfd_decomposition(x, 'tvfilt', 2);
```

where `y_comps` is a cell array of the 2 components and `y` is sum of the components.

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
which sets the TFD to be of size `256 x 32,768`. Accuracy for the xTFD method is tied to
the `Nfreq`: a larger value reduces error but increase computational load and memory
requirements. The default value for `Nfreq` is set to 8,192.


## 2-tone test
To run the 2-tone test for each method, call:

```matlab
time_mask = false;

sep_2tones_test('xtfd', time_mask);

```

where `time_mask` is set to `true` to use the time-separated signal components. Need to
call `sep_2tones_test` for each method. That is, for all methods:

```matlab
time_mask = false;

all_methods = {'tvfilt', 'xtfd', 'vmd', 'ssst', 'edf', 'tvemd'};
     
for n = 1:length(all_methods)
	sep_2tones_test(all_methods{n}, time_mask, true, false, 50);
end

% and plot all:
plot_2tone_tests_compare(false);
```

and set `time_mask = true;` and same again.

## Test signals

### a) bat signal

To generate the plots for bat signals, decomposing to a maximum of 5 components, do as follows:
```matlab
     print_ = false;
     all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'msst', 'vmd'};
     
     for n = 1:length(all_methods)
          compare_methods_generic_signal_plot(all_methods{n}, 'bat', print_);
     end
```


### b) fractional Gaussian noise

```matlab
       print_ = false;
       all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'msst', 'vmd'};

       for n = 1:length(all_methods)
           compare_methods_noise_signals_plot(all_methods{n}, print_);
       end
```

### c) multicomponent nonlinear FM signal with noise
```matlab
	compare_methods_generic_signal_plot('tvfilt', 'nnlfm4', false);
	compare_methods_generic_signal_plot('xtfd', 'nnlfm4', false);
```

and plot:

```matlab
	plot_nnlfm4_testsignal('nnlfm4', false);
```
