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
