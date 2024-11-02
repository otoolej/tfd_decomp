function plot_bat_tfd(print_)
    if(nargin < 1) print_ = false; end

    [x, x_components, Fs, all_params] = set_signal_parameters('bat', false);
    N = length(x);
    qtfd = qtfd_sep_kern(x, {63, 'hann'}, {63, 'dolph', 100}, N, N);
    image_plot(qtfd, 1, 'Signal');
    if(print_)
        print2eps('pics/bat_test/bat_TFD.eps')
    end
    

    [x, x_components, y, y_comps] = compare_methods_testsignals('bat', {'xtfd'}, false);
    qtfd = zeros(N, N);
    for n = 1:4
        qtfd_tmp = qtfd_sep_kern(y_comps{n}, {63, 'hann'}, {63, 'dolph', 100}, N, N);
        qtfd = qtfd + qtfd_tmp;
    end
    image_plot(qtfd, 2, 'XTFD');
    if(print_)
        print2eps('pics/bat_test/bat_TFD_xtfd_sum_all_components.eps')
    end
    

    [x, x_components, y, y_comps] = compare_methods_testsignals('bat', {'tvfilt'}, false);
    qtfd = zeros(N, N);
    for n = 1:4
        qtfd_tmp = qtfd_sep_kern(y_comps{n}, {63, 'hann'}, {63, 'dolph', 100}, N, N);
        qtfd = qtfd + qtfd_tmp;
    end
    image_plot(qtfd, 3, 'TV-FILT');
    if(print_)
        print2eps('pics/bat_test/bat_TFD_tvfilt_sum_all_components.eps')
    end

end


function image_plot(qtfd, fig_num, title_str)
    FONT_NAME = 'helvetica';
    FONT_SIZE = 22;

    qtfd = qtfd ./ max(qtfd(:));
    qtfd(qtfd(:) > 0.3) = 0.3;


    set_figure(fig_num);
    pp = get(gcf, 'position');
    set(gcf, 'position', [pp(1:2) 570 400]);
    cc = cubehelix(256, 0.5, -1.5, 1, 1, [0.0, 1.0]);
    cc = flipud(cc);
    colormap(cc); % , 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]));

    imagesc(qtfd);
    axis('xy'); axis('tight'); axis('square');
    hx = gca();
    % hx.XAxis.Visible = 'off';
    % hx.YAxis.Visible = 'off';
    hx.Box = 'on';
    hx.XAxis.TickValues = [];
    hx.YAxis.TickValues = [];
    xlabel('frequency', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
    ylabel('time', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
    title(title_str, 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
    % set_gca_fonts(FONT_NAME, FONT_SIZE, hx);

end
