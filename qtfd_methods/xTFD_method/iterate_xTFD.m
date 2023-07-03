%-------------------------------------------------------------------------------
% iterate_xTFD: iterate the xTFD method
%
% Syntax: [y, x_residual, y_components] = iterate_xTFD(x, Fs, params, N_components, N_iters)
%
% Inputs: 
%     x, Fs, params, N_components, N_iters - 
%
% Outputs: 
%     [y, x_residual, y_components] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 06-09-2021
%
% last update: Time-stamp: <2021-09-21 11:58:21 (otoolej)>
%-------------------------------------------------------------------------------
function [y, y_components] = iterate_xTFD(x, Fs, params, N_components, N_iters, db_plot)
if(nargin < 2 || isempty(Fs)), Fs = 1; end
if(nargin < 3 || isempty(params)), params = decomp_params; end
if(nargin < 4 || isempty(N_components)), N_components = 1; end
if(nargin < 5 || isempty(N_iters)), N_iters = N_components; end
if(nargin < 6 || isempty(db_plot)), db_plot = false; end



x_residual = x;
for n = 1:N_iters

    [y, x_res] = extract_components_xTFD(x_residual, Fs, params, N_components, db_plot);

    if(db_plot)
        disp('--- paused; hit key to continue ---'); pause;
    end


    y_components{n} = y;
    if(~isempty(y))
        x_residual = x_residual - y;         
    end
end

y = nansum( cat(1, y_components{:}), 1);


% db_plot = false;
% if(db_plot)
%     plot_components(x, y_components, Fs, 29, N_components, false);
% end
