%-------------------------------------------------------------------------------
% compare_methods_testsignals: 
%
% Syntax: [] = compare_methods_testsignals(signal_name, methods_subset)
%
% Inputs: 
%     signal_name, methods_subset - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 03-09-2021
%
% last update: Time-stamp: <2023-07-06 07:45:10 (otoolej)>
%-------------------------------------------------------------------------------
function [x, x_components, y, y_comps] = compare_methods_testsignals(signal_name, methods_subset, dbplot)
if(nargin < 1 || isempty(signal_name)), signal_name = 'lfm1'; end
if(nargin < 2 || isempty(methods_subset)), methods_subset = {'xTFD', 'WSST'}; end
if(nargin < 3 || isempty(dbplot)), dbplot = false; end



%---------------------------------------------------------------------
% 1. load signals
%---------------------------------------------------------------------
[x_nlfm, x_nlfm_comps, x_lfm_comps] = nlfm_test_signals(false);


%---------------------------------------------------------------------
% 2. default parameters for all methods
%---------------------------------------------------------------------
[x, x_components, Fs, all_params] = set_signal_parameters(signal_name, false);

N_components = length(x_components);

% need to keep original sampling rates:
UPSAMPLE = false;
if(UPSAMPLE)
    xo = x;
    xco = x_components;
end



%---------------------------------------------------------------------
% 2. compare methods
%---------------------------------------------------------------------
db_plot = true;
mse_all = [];
for n = 1:length(methods_subset)
    % set_figure(8 + n);
    % plot(y_comps{1} - x_components{1}'); % plot(y_comps{2} - x_components{2}');

    % upsampling for specific method:
    if(UPSAMPLE)
        switch methods_subset{n}
          case {'ITD', 'EMD'}
            x = resample(x, 5, 1);
            for p = 1:length(x_components)
                x_components{p} = resample(x_components{p}, 5, 1);
            end
        end
    end

    % make sure components are ordered correctly:
    p_st = all_params(find(strcmp({all_params.method}, lower(methods_subset{n}))));
    if(isempty(p_st))
        params = [];
    else
        params = p_st.params;
    end
    
    % dispVars(methods_subset{n});
    % disp(p_st.params);

    [y, y_comps] = decomp_all_methods(x, Fs, methods_subset{n}, N_components, params, dbplot);

    y = y(:);
    if(db_plot && ~isempty(y))
        err_plot(x, y_comps, y, 2000 + n * 100, methods_subset{n});
    end

    
    mse_tmp = [];
    rerr_tmp = [];
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
        rerr_tmp = [rerr_tmp (diff_l2 / comp_l2)];                
        rr = corrcoef(x_components{isort(1)}(:), yc(:));
        mse_tmp = [mse_tmp rr(1, 2)];

        
    end
    % dispVars(mse_tmp);
    mse_all{n} = mean(mse_tmp);
    re_all{n} = mean(rerr_tmp);
    snr_all{n} = 10 * log10(mean(1 ./ rerr_tmp));

    % restore original signal
    if(UPSAMPLE)
        switch methods_subset{n}
          case {'ITD', 'EMD'}
            x = xo;
            x_components = xco;
        end
    end
end

fprintf('signal %s:\n', signal_name);
print_table([mse_all{:}; re_all{:}; snr_all{:}]', {'CORR', 'RE', 'SNR'}, ...
            methods_subset, [], [], [4]);
% print_table(mse_all', {'MSE'}, methods_subset);



