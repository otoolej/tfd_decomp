# MATLAB golden fixture

Run `tests/matlab/generate_reference_fixtures.m` from MATLAB after adding this repository
to the path. It writes `matlab_reference.mat` and `matlab_dataset_reference.mat` here. The
first contains staged intermediate results for a synthetic chirp; the second contains final
outputs for compact segments of all bundled signal types. Python golden-parity tests skip
when the generated files are absent.
