%-------------------------------------------------------------------------------
% iterative_tf_filt_JOT: 
%
% Syntax: y = iterative_tf_filt_JOT(x, params, N_components, N_iter, db_plot)
%
% Inputs: 
%     x, params, N_components, N_iter, db_plot - 
%
% Outputs: 
%     y - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 21-09-2021
%
% last update: Time-stamp: <2021-09-21 12:03:44 (otoolej)>
%-------------------------------------------------------------------------------
function [y, y_components] = iterate_tv_filt_JOT(x, params, N_components, N_iter, db_plot)
if(nargin < 2 || isempty(params)), params = tvfilt_params(length(x)); end
if(nargin < 4 || isempty(N_iter)), N_iter = N_components; end
if(nargin < 5 || isempty(db_plot)), db_plot = false; end




x_residual = x;
for n = 1:N_iter

    y = tf_decomposition_v1_JOT(x_residual, params, N_components, db_plot);

    if(db_plot)
        disp('--- paused; hit key to continue ---'); pause;
    end


    y_components{n} = y;
    if(~isempty(y))
        x_residual = x_residual - y;         
    end
end

y = nansum( cat(1, y_components{:}), 1);

