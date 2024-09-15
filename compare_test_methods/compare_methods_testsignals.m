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
% last update: Time-stamp: <2024-09-15 20:51:30 (otoolej)>
%-------------------------------------------------------------------------------
function [x, x_components, y, y_comps, comp_time] = compare_methods_testsignals(signal_name, methods_subset, dbplot)
if(nargin < 1 || isempty(signal_name)), signal_name = 'lfm1'; end
if(nargin < 2 || isempty(methods_subset)), methods_subset = {'xTFD', 'WSST'}; end
if(nargin < 3 || isempty(dbplot)), dbplot = false; end



%---------------------------------------------------------------------
% 1. load the the signals and default parameters for all methods
%---------------------------------------------------------------------
[x, x_components, Fs, all_params] = set_signal_parameters(signal_name, false);

N_components = length(x_components);

%  for testing how many signal components with <90% of signal:
N_components = 100;

% need to keep original sampling rates:
UPSAMPLE = false;
if(UPSAMPLE)
    xo = x;
    xco = x_components;
end



%---------------------------------------------------------------------
% 2. compare methods
%---------------------------------------------------------------------
db_plot = false;
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
    
    dispVars(methods_subset{n});
    dispVars(N_components);
    % disp(p_st.params);
    
    [y, y_comps, comp_time] = decomp_all_methods(...
        x, Fs, methods_subset{n}, N_components, params, dbplot);

    y = y(:);
    if(db_plot && ~isempty(y))
        err_plot(x, y_comps, y, 2000 + n * 100, methods_subset{n});
    end

    switch signal_name
      case {'bat', 'noise', 'white-noise'}
        re_all{n} = sqrt(sum(abs(x - y) .^ 2)) / sqrt(sum(abs(x) .^ 2));
        snr_all{n} = 10 * log10(mean(1 ./ re_all{n}));
        r_tmp = corrcoef(x, y);
        mse_all{n} = r_tmp(1, 2);
        
      otherwise

        mse_tmp = [];
        rerr_tmp = [];
        rerr_tmp = [];
        if(strcmp(signal_name, 'nnlfm4'))
            l_end = 3;
        else
            l_end = length(y_comps);
        end
        for p = 1:l_end
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
            if(var(yc) == 0)
                rr(1, 2) = 0;
            end
            mse_tmp = [mse_tmp rr(1, 2)];
            dispVars(mse_tmp);
            % set_figure(30); 
            % plot(x_components{isort(1)});            
            % plot(yc);
            % input('any key..')
        end
        % dispVars(mse_tmp);
        mse_all{n} = mean(mse_tmp);
        re_all{n} = mean(rerr_tmp);
        snr_all{n} = 10 * log10(mean(1 ./ rerr_tmp));
    end
    
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
            methods_subset, [], [], [3]);
% print_table(mse_all', {'MSE'}, methods_subset);



