# Signal Decomposition using Quadratic Time--Frequency Distributions (QTFDs)

Two methods of signal decomposition using QTFDs:
1. time-varying filtering approach
2. sinusoidal model approach with cross-TFD phase correction

## directory structure

```
.
├── data
│   └── test_signals    # signals for testing/comparing the methods
├── docs
│   └── papers          # relevant papers
├── existing_methods    # code for existing methods (e.g. EMD)
│   ├── emd
│   ├── tv_emd
│   ├── vmd
│   └── vncmd
├── extract_IFs         # common methods for extracting the IFs and IAs
├── README.md
├── tvfilt_decomp_NS    # METHOD-1
└── xtfd_decomp_JOT     # METHOD-2
```

# THINGS NATHAN HAS DONE: 29/08/19
Uploaded first go at algorithm: file of interest is
`tvfilt_decomp_NS\tf_decomposition_demo_v1.m`, the remaining files are called by this demo
Uploaded EEG test signal into data\test_signals: `file test_signal_eeg1.mat`

## Comments: 
in general pretty good but weird missing bit in reconstructed signal at around 1750 samples
(edge linking not picking up the LF component at this point in time). I will look into
it. Good things to note is that it shows temporal segmentation, not just frequency based
segmentation.

# THINGS NATHAN HAS DONE: 12/07/21
I have improved the edge-linking algorithm (it now tracks bifurcations and selects the best version at a time). Also mucking about with filter bandwidths in the decomposition step (thinking of setting two options - full band or narrowband decomposition, where full band decomposes the whole TF plane so reconstruction will be approximately perfect, if that makes any sense, and narrowband where only the main components are extracted). I have implemented the process on the three test signals and get something else - although it is pretty dependent on the parameter values set within the subfunctions. I wonder if there are some decent defaults here that might makes sense (sqrt(N) is always a good one). 


