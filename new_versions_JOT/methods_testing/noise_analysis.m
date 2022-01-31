%-------------------------------------------------------------------------------
% noise_analysis: 
%
% Syntax: n_st = noise_analysis(signal_name, methods_subset)
%
% Inputs: 
%     signal_name, methods_subset - 
%
% Outputs: 
%     n_st - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 24-09-2021
%
% last update: Time-stamp: <2021-09-26 20:57:52 (otoolej)>
%-------------------------------------------------------------------------------
function n_st = noise_analysis(signal_name, methods_subset, db_plot)
if(nargin < 1 || isempty(signal_name)), signal_name = 'lfm1'; end
if(nargin < 2 || isempty(methods_subset)), methods_subset = {'xTFD', 'WSST'}; end
if(nargin < 3 || isempty(db_plot)), db_plot = false; end


%---------------------------------------------------------------------
% 0. set the noise parameters
%---------------------------------------------------------------------
noise_type = 'wgn';
snr_ranges = -3:15;
L_snrs = length(snr_ranges);


%---------------------------------------------------------------------
% 1. load signals (put may need more?)
%---------------------------------------------------------------------
[x, x_components, Fs, all_params] = select_signal_withparams(signal_name, db_plot);
N = length(x);
N_components = length(x_components);

switch lower(noise_type)
  case 'wgn'
    w = zscore(randn(N, 1));
end


re_snrs = zeros(L_snrs, length(methods_subset));

for n = 1:L_snrs
    xw = x + w ./ 10^(snr_ranges(n)/10);
    dispVars(snr_ranges(n), 10^(snr_ranges(n)/10))

    % mse_all = []; rerr_all = []; snr_all = [];
    for m = 1:length(methods_subset)
        % make sure components are ordered correctly:
        p_st = all_params(find(strcmp({all_params.method}, lower(methods_subset{m}))));
        if(isempty(p_st))
            params = [];
        else
            params = p_st.params;
        end
        
        % dispVars(methods_subset{n});
        % disp(p_st.params);

        [y, y_comps] = decomp_all_methods(xw, Fs, methods_subset{m}, N_components, params, db_plot);
        % dispVars(size(y), size(y_comps));

        if(size(y, 2) > 1)
            y = y(:);
        end
        if(db_plot && ~isempty(y_comps))
            err_plot(xw, y_comps, y, 2000 + (m * 100), methods_subset{m});
        end

        [mse_all{m}, re_all{m}, snr_all{m}] = cal_err(x_components, y_comps);

        re_snrs(n, m) = re_all{m};
    end
    % fprintf('signal %s:\n', signal_name);
    print_table([mse_all{:}; re_all{:}; snr_all{:}]', {'MSE', 'RE', 'SNR'}, ...
                methods_subset, [], [], [4]);
    % if(n == L_snrs), keyboard; end
    
end

fprintf('\n');
print_table([snr_ranges' re_snrs], ['SNR' methods_subset], [], [], [], ...
            [1 repmat(3, 1, length(methods_subset))]);
