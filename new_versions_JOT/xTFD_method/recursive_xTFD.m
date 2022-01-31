%-------------------------------------------------------------------------------
% recursive_xTFD: iterate the xTFD method
%
% Syntax: [y, x_residual, x_components] = recursive_xTFD(x, Fs, params, N_components, N_iters)
%
% Inputs: 
%     x, Fs, params, N_components, N_iters - 
%
% Outputs: 
%     [y, x_residual, x_components] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 06-09-2021
%
% last update: Time-stamp: <2021-09-06 17:37:38 (otoolej)>
%-------------------------------------------------------------------------------
function [y, x_residual, x_components] = recursive_xTFD(x, Fs, params, N_components, N_iters)
if(nargin < 2 || isempty(Fs)), Fs = 1; end
if(nargin < 3 || isempty(params)), params = decomp_params; end
if(nargin < 4 || isempty(N_components)), N_components = 1; end
if(nargin < 5 || isempty(N_iters)), N_iters = N_components; end


x_residual = x;
for n = 1:N_iters

    [y, x_res] = extract_components_xTFD(x_residual, Fs, params, 1, true);
    disp('--- paused; hit key to continue ---'); pause;


    x_components{n} = y;
    if(~isempty(y))
        x_residual = x_residual - y;         
    end
end



db_plot = true;
if(db_plot)
    plot_components(x, x_components, Fs, 29, N_components, false);
end
