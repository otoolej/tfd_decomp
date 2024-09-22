all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'vmd', 'vncmd'};
[x, x_components, y, y_comps] = compare_methods_testsignals('bat', all_methods, false);
[x, x_components, y, y_comps] = compare_methods_testsignals('nnlfm4', all_methods, false);
[x, x_components, y, y_comps] = compare_methods_testsignals('noise', all_methods, false);




[x, x_components, y, y_comps] = compare_methods_testsignals('bat', {'tvfilt'}, true);
