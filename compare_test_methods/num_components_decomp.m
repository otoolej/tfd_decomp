%-------------------------------------------------------------------------------
% num_components_decomp: 
%
% Syntax: [] = num_components_decomp()
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
% Started: 09-07-2023
%
% last update: Time-stamp: <2023-07-09 13:42:13 (otoolej)>
%-------------------------------------------------------------------------------
function [] = num_components_decomp()

    
all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'vmd', 'vncmd'};
signal_types = {'2tone1', '2tone2', 'nnlfm4', 'noise', 'bat'};
% signal_types = {'2tone2'};
% all_methods = {'vncmd'};

num_comps = zeros(length(signal_types), length(all_methods));

for p = 1:length(signal_types)
    for n = 1:length(all_methods)
        [x, x_components, y, y_comps] = compare_methods_testsignals(...
            signal_types{p}, all_methods(n), false);
        
        x_eng = sum(abs(x) .^ 2);
        q = 1;
        y_eng = sum(abs(y_comps{1}) .^ 2);
        if(y_eng <  (0.9 * x_eng))
            q = 1;
            while (y_eng < (0.9 * x_eng) && ~isnan(y_eng))
                q = q + 1;
                if(q > length(y_comps))
                    q = q - 1;
                    break;
                end
                y_eng = sum(abs(sum(cat(1, y_comps{1:q}))) .^ 2);
                dispVars(x_eng, y_eng, q);
            end
        end
        if(y_eng < (0.9 * x_eng))
            q = 0;
        end
        num_comps(p, n) = q; 
    end
end

print_table(num_comps, all_methods, signal_types, [6], [], 0);


