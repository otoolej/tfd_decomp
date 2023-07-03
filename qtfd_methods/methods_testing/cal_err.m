%-------------------------------------------------------------------------------
% cal_err: 
%
% Syntax: err_st = cal_err(x_components, y_comps)
%
% Inputs: 
%     x_components, y_comps - 
%
% Outputs: 
%     err_st - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 24-09-2021
%
% last update: Time-stamp: <2021-09-24 17:56:06 (otoolej)>
%-------------------------------------------------------------------------------
function [mse_all, re_all, snr_all] = cal_err(x_components, y_comps, display)
if(nargin < 3 || isempty(display)), display = true; end


mse_tmp = [];
rerr_tmp = [];
for p = 1:length(y_comps)
    yc = y_comps{p};
    yc(isnan(yc)) = 0;

    r = zeros(1, length(x_components));
    for q = 1:length(x_components)
        r(q) = corr(x_components{q}(:), yc(:));
    end
    [~, isort] = max(r);

    diff_l2 = sqrt(sum(abs(x_components{isort(1)}(:)' - yc(:).') .^ 2));
    comp_l2 = sqrt(sum(abs(x_components{isort(1)}) .^ 2));
    rr = corrcoef(x_components{isort(1)}(:), yc(:));
    mse_tmp = [mse_tmp rr(1, 2)];
    rerr_tmp = [rerr_tmp (diff_l2 / comp_l2)];        
end

% dispVars(mse_tmp);
mse_all = mean(mse_tmp);
re_all = mean(rerr_tmp);
snr_all = 10 * log10(mean(1 ./ rerr_tmp));        
