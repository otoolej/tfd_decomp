%-------------------------------------------------------------------------------
% tfd_decomposition: wrapper function for TV-filt and xTFD decomposition methods
%
% Syntax: [y, y_comps] = tfd_decomposition(x, method, N_components, params, db_plot)
%
% Inputs: 
%     x            - signal (length-N)
%     method       - decomposition method; either 'xtfd' or 'tvfilt' (default)
%     N_components - maximum number of components to extract (default=10)
%     params       - parameters for the methods (object; see 'decomp_params.m')
%     db_plot      - plot (default=false)
% 
%
% Outputs: 
%     y       - output signal consisting of sum of components
%     y_comps - cell array containing the components
%
%

% John M. O' Toole, University College Cork
% Started: 20-04-2022
%
% last update: Time-stamp: <2025-06-12 21:02:36 (otoolej)>
%-------------------------------------------------------------------------------
function [y, y_comps] = tfd_decomposition(x, method, N_components, params, db_plot)
if(nargin < 2 || isempty(method)), method = 'tvfilt'; end
if(nargin < 3 || isempty(N_components)), N_components = 10; end
if(nargin < 4 || isempty(params))
    params = decomp_params(length(x), method);
end
if(nargin < 5 || isempty(db_plot)), db_plot = false; end


% must be a row vector and even length signal:
x = x(:)';
x = x(1:end - mod(length(x), 2));



% if preprocessing the signal:
if(params.pad_signal)
    N = length(x);
    No = N;
    xo = x;
    x = preprocessing_signal(x, true);
    Lpad = length(x) - N;
    N = length(x);
end



switch upper(method)
  case {'XTFD'}
    %---------------------------------------------------------------------
    % cross-TFD method
    %---------------------------------------------------------------------
    [y, ~, y_comps] = extract_components_xTFD(x, 1, params, N_components, db_plot);


  case {'TVFILT'}
    %---------------------------------------------------------------------
    % time-varying filtering method
    %---------------------------------------------------------------------
    d = tf_decomposition_tvfilt(x, params, N_components, db_plot);

    if(size(d, 1) < N_components)
        N_components = size(d, 1);
    end
    d = d(1:N_components, :);

    y_comps = num2cell(d, 2);
    y = nansum(d, 1).';


  otherwise
    erro('Unknown method: should be either tvfilt or xtfd');
end





if(params.pad_signal)
    y = y(floor(Lpad / 2) + 1:end - floor(Lpad / 2));
    for p = 1:length(y_comps)
        y_comps{p} = y_comps{p}(floor(Lpad / 2) + 1:end - floor(Lpad / 2));
    end

    x = xo;
end

if(db_plot && ~isempty(y_comps))
    err_plot(y, y_comps, x, 20, method);
end
