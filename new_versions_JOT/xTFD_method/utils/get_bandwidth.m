%-------------------------------------------------------------------------------
% get_bandwidth: estimate bandwidth of component 
%
% Syntax: [] = get_bandwidth(S_tslice, n_peak, amp_thres, db_plot)
%
% Inputs: 
%     S_tslice  - time-slice of the TFD
%     n_peak    - index of peak
%     amp_thres - amplitude value to measure bandwidth (default = peak/2)
%     db_plot   - plot (bool)
%
% Outputs: 
%     bw - bandwidth 
%


% John M. O' Toole, University College Cork
% Started: 13-07-2021
%
% last update: Time-stamp: <2021-09-02 15:20:41 (otoolej)>
%-------------------------------------------------------------------------------
function [bw] = get_bandwidth(S_tslice, n_peak, amp_thres, dbplot)
if(nargin < 2 || isempty(n_peak)), n_peak = []; end
if(nargin < 3 || isempty(amp_thres)), amp_thres = 0.5; end
if(nargin < 4 || isempty(dbplot)), dbplot = false; end

if(isempty(n_peak))
    [~, n_peak] = max(S_tslice);
end

% maximum point:
amp_x = S_tslice(n_peak) * amp_thres;

% find closest points to the mark:
ileft = find(S_tslice(1:n_peak) <= amp_x, 1, 'last');
iright = n_peak + find(S_tslice((n_peak + 1):end) <= amp_x, 1, 'first');

if(isempty(ileft))
    ileft = 1;
end
if(isempty(iright))
    iright = length(S_tslice);
end

%---------------------------------------------------------------------
% linear interpolation of estimate
%---------------------------------------------------------------------
if(S_tslice(ileft) < amp_x)
    [m, c] = slope_intercep_line(S_tslice, ileft, 1);

    ileft = (amp_x - c) / m;
end
if(S_tslice(iright) < amp_x)
    [m, c] = slope_intercep_line(S_tslice, iright, -1);

    iright = (amp_x - c) / m;
end

% bandwidth:
bw = iright - ileft;

if(bw < 0)
    bw = [];
end


if(dbplot)
    set_figure(67);
    plot(S_tslice_interp, '-o');
    ys = ylim;
    xs = xlim;
    line([1, 1] .* ileft, ys, 'color', 'k');
    line([1, 1] .* iright, ys, 'color', 'k');    
    line(xs, [1, 1] .* amp_x, 'color', 'r');
end




function [m, c] = slope_intercep_line(S_tslice, idx, ishift)
%---------------------------------------------------------------------
% get slope and intercept of line
%---------------------------------------------------------------------
m = (S_tslice(idx + ishift) - S_tslice(idx)) / ishift;
c = S_tslice(idx) - (m * idx);



