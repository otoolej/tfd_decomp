%-------------------------------------------------------------------------------
% num_components_decomp: minimum number of combined components to reach the 90% 
% energy threshold
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
% last update: Time-stamp: <2024-10-29 23:55:09 (otoolej)>
%-------------------------------------------------------------------------------
function [num_comps, all_methods, signal_types] = num_components_decomp(all_methods, signal_types)
if(nargin < 1 || isempty(all_methods))
    all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'vmd', 'vncmd'};
end
if(nargin < 2 || isempty(signal_types))
    signal_types = {'2tone1', '2tone2', 'nnlfm4', 'white-noise','noise', 'bat'};
end



num_comps = zeros(length(signal_types), length(all_methods));
corr_x = zeros(length(signal_types), length(all_methods));

for p = 1:length(signal_types)
    for n = 1:length(all_methods)
        [x, x_components, y, y_comps, ~] = compare_methods_testsignals(...
            signal_types{p}, all_methods(n), false);
        dispVars(size(y_comps));
        
        x_eng = sum(abs(x) .^ 2);
        q = 1;
        y_est = y_comps{1}(:)';        
        y_eng = sum(abs(y_est) .^ 2);
        if(y_eng <  (0.9 * x_eng))
            q = 1;
            while (y_eng < (0.9 * x_eng) && ~isnan(y_eng))
                q = q + 1;
                if(q > length(y_comps))
                    q = q - 1;
                    break;
                end
                y_est = sum_components(y_comps(1:q));
                y_eng = sum(abs(y_est) .^ 2);
                
                dispVars(q, x_eng, y_eng, 100 * y_eng / x_eng);
                
                % dispVars(size(y_est'), size(x));
                % corr_x(p, n) = corr(y_est', x);
                
                % dispVars(x_eng, y_eng, q);
            end
        end
        if(y_eng < (0.9 * x_eng))
            q = -q;
        end
        dispVars(y_eng, x_eng, q, 100 * y_eng / x_eng);
        
        corr_x(p, n) = corr(y_est', x);        
        num_comps(p, n) = q; 
    end
end

print_table(num_comps, all_methods, signal_types, [6], [], 0);
print_table(corr_x, all_methods, signal_types, [6], [], 2);


function y = sum_components(y_comps)
%---------------------------------------------------------------------
% sum components
%---------------------------------------------------------------------
if(iscell(y_comps))
    y = 0;
    for n = 1:length(y_comps)
        y = y + y_comps{n}(:)';
    end
else
    y = y_comps;
end


