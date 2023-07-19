%-------------------------------------------------------------------------------
% timeit_methods: computation time for methods
%
% Syntax: [] = timeit_methods()
%
% Inputs: 
%      - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 19-07-2023
%
% last update: Time-stamp: <2023-07-19 19:46:59 (otoolej)>
%-------------------------------------------------------------------------------
function [comp_time_mean, comp_time_std] = timeit_methods()

all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'vmd', 'vncmd'};
signal_types = {'2tone1', '2tone2', 'nnlfm4', 'noise', 'bat'};
% signal_types = {'2tone2'};
% all_methods = {'vncmd'};

% all_methods = {'tvfilt'};
n_iter = 100;

comp_time_mean = zeros(length(signal_types), length(all_methods));
comp_time_std = zeros(length(signal_types), length(all_methods));

for p = 1:length(signal_types)
    for n = 1:length(all_methods)
        
        cp = zeros(1, n_iter);
        for m = 1:n_iter
            [~, ~, ~, ~, cp(m)] = compare_methods_testsignals(...
                signal_types{p}, all_methods(n), false);
        end
            
        comp_time_mean(p, n) = mean(cp);
        comp_time_std(p, n) = std(cp);        
    end
end

print_table(comp_time_mean, all_methods, signal_types, [8], [], 4);
print_table(comp_time_std, all_methods, signal_types, [10], [], 6);
