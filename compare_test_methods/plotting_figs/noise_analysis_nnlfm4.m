%-------------------------------------------------------------------------------
% noise_analysis_nnlfm4: 
%
% Syntax:  = noise_analysis_nnlfm4(all_methods)
%
% Inputs: 
%     all_methods - 
%
% Outputs: 
%      - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 09-09-2024
%
% last update: Time-stamp: <2024-09-15 22:00:46 (otoolej)>
%-------------------------------------------------------------------------------
function all_out_tb = noise_analysis_nnlfm4(plot_only)
    if(nargin < 1) plot_only = false; end
    all_methods = {'xtfd', 'tvfilt', 'efd', 'tvemd', 'ssst', 'vmd', 'vncmd'};
    % all_methods = {'ssst'};
    
    if(plot_only)
        d = load('data/plots/noise_analysis.mat');
        all_out_tb = d.out;
    else

        all_out_tb = table();
        for n = 1:length(all_methods)
            out_tmp = noise_analysis_per_method(all_methods{n});
            all_out_tb = [all_out_tb; out_tmp];
        end
    end
    plot_snr_all_methods(all_out_tb, all_methods, 'corr',1);
    plot_snr_all_methods(all_out_tb, all_methods, 're', 2);
end



function out_tb = noise_analysis_per_method(method_str)
    if(nargin < 1) method_str = 'tvfilt'; end


    [x, x_components, Fs, all_params] = set_signal_parameters('nn-nnlfm4', false);

    d = load('data/test_signals/ffgn_1_02_1_512_0_signal.mat');
    n = d.x(1:256);
    N_components = length(x_components);
    dispVars(N_components);

    db_range = [-4, -2, 0, 2, 4, 6, 8, 10];

    p_st = all_params(find(strcmp({all_params.method}, lower(method_str))));
    if(isempty(p_st))
        params = [];
    else
        params = p_st.params;
    end

    n_db = length(db_range);
    mse_all = zeros(1, n_db); re_all = zeros(1, n_db); snr_all = zeros(1, n_db);
    q = 1;
    db_plot = false;
    for dc = db_range
        n_scalar = sqrt((sum(abs(x) .^ 2) / sum(abs(n) .^ 2)) * 10^(-dc/10));
        ns = n .* n_scalar;
        fprintf('SNR=%g (dB)\n', 10 * log10(sum(abs(x) .^ 2) / sum(abs(ns) .^ 2)));
        x_test = x + ns(:);

        
        [y, y_comps, comp_time] = decomp_all_methods(...
            x_test, Fs, method_str, N_components, params, db_plot);
        y = y(:);
        if(db_plot)
            input('any key. ..');
        end
        l_end = 3;
        mse_tmp = zeros(l_end, 1);
        rerr_tmp = zeros(l_end, 1);
        for p = 1:l_end
            yc = y_comps{p};
            yc(isnan(yc)) = 0;

            % order the estimated components to the highest correlation with real components:
            r = zeros(1, length(x_components));
            for m = 1:length(x_components)
                r(m) = corr(x_components{m}(:), yc(:));
            end
            [~, isort] = max(r);

            diff_l2 = sqrt(sum(abs(x_components{isort(1)}(:)' - yc(:).') .^ 2));
            comp_l2 = sqrt(sum(abs(x_components{isort(1)}) .^ 2));
            rerr_tmp(p) = (diff_l2 / comp_l2);                
            rr = corrcoef(x_components{isort(1)}(:), yc(:));
            if(var(yc) == 0)
                rr(1, 2) = 0;
            end
            mse_tmp(p) = rr(1, 2);
            % dispVars(mse_tmp);
            % set_figure(30); 
            % plot(x_components{isort(1)});            
            % plot(yc);
            % input('any key..');
        end
        mse_all(q) = mean(mse_tmp);
        re_all(q) = mean(rerr_tmp);
        snr_all(q) = 10 * log10(mean(1 ./ rerr_tmp));
        db_str{q} = sprintf('%d dB', dc);
        q = q + 1;
    end
    print_table([mse_all; re_all; snr_all]', {'CORR', 'RE', 'SNR'}, db_str, [], [], [3]);

    out = [mse_all; re_all; snr_all]';
    out_tb = cell2table([num2cell(db_range)' num2cell(out)], 'Variablenames', {'SNR', 'corr', 'RE', 'snr_re'});
    out_tb.method = repmat({method_str}, height(out_tb), 1);
end


function str_names = convert_str_method_names(all_methods)
    str_names = {};
    for n = 1:length(all_methods)
        switch all_methods{n}
          case {'xtfd', 'efd', 'vmd'}
            fmt_str = upper(all_methods{n});
          case 'tvfilt'
            fmt_str = 'TV-FILT';
          case 'tvemd'
            fmt_str = 'TV-EMD';
          case 'ssst'
            fmt_str = 'SSFT';
          case 'vncmd'
            fmt_str = 'VN-CMD';
          otherwise
            error('??');
        end
        str_names{n} = fmt_str;
    end
end



function plot_snr_all_methods(out_tb, methods, err_type, fignum)
% Extract unique methods and markers
    markers = {'o', 's', 'd', '^', 'hexagram', 'v', '<', '>'}; % Different markers for each method
    lc = cubehelix(8, 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]);
    llcube = cubehelix(10);
    % llcube = llcube([1, 7, 2, 6, 3, 5, 4], :);
    llcube = llcube([2, 3:8], :);
    % methods = fliplr(methods);
    str_names = convert_str_method_names(methods);
    FONT_NAME = 'helvetica';
    FONT_SIZE = 12;
    

    % Create a new figure
    fig = set_figure(fignum);
    xy = fig.Position;
    set(fig, 'position', [xy(1), xy(2), 570, 366]);

    % Loop through each method and plot the corresponding data
    for i = 1:length(methods)
        idx = strcmp(out_tb.method, methods{i}); % Get indices for the current method
        if(strcmp(err_type, 'corr'))
            y = out_tb.corr(idx);
        else
            y = out_tb.RE(idx);
        end

        plot(out_tb.SNR(idx), y, 'DisplayName', str_names{i}, ...
             'Marker', markers{i}, 'LineWidth', 0.8, 'color', llcube(i, :), ...
             'markerfacecolor', llcube(i, :), 'MarkerSize', 5, ...
             'markeredgecolor', llcube(i, :));
    end
    set_gca_fonts(FONT_NAME, FONT_SIZE - 2, gca());
    % Add labels and legend
    xlabel('signal to noise ratio (dB)', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
    ylabel('correlation', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
    leg = legend('Position', [0.7310, 0.2556, 0.1554, 0.2167]);

    if(strcmp(err_type, 'corr'))
        leg = legend('Position', [0.5028 0.2445 0.1768 0.2486]);
        ylabel('correlation', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
        fname = 'pics/noise_analysis_corr_v1.eps';
    else
        ylabel('relative error', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
        leg = legend('position', [0.1464, 0.1276, 0.1737, 0.2486]);
        fname = 'pics/noise_analysis_re_v1.eps';
    end
    print2eps(fname);
end
