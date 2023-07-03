%-------------------------------------------------------------------------------
% compare_methods_noise_signals_plot: compare a few methods using fractional noise
%
% Syntax: [] = compare_methods_noise_signals_plot()
%
% Inputs: 
%      - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 12-11-2021
%
% last update: Time-stamp: <2023-07-02 09:51:23 (otoolej)>
%-------------------------------------------------------------------------------
function compare_methods_generic_signal_plot(method_str, signal_type, print_)
if(nargin < 1 || isempty(method_str)), method_str = 'vmd'; end
if(nargin < 2 || isempty(signal_type)), signal_type = 'bat'; end
if(nargin < 3 || isempty(print_)), print_ = false; end


scale_fig = true;


FONT_NAME = 'helvetica';
FONT_SIZE = 14;
if(scale_fig)
    FONT_SIZE = 14;
end


all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'msst', 'vmd'};

[x, x_components, y, y_comps] = compare_methods_testsignals(signal_type, {method_str}, false);

%---------------------------------------------------------------------
% 
%---------------------------------------------------------------------
 
hp = plot_components(x, y_comps, 1, 900, length(y_comps), false);
% if(length(hp) > 20)
%     yticklabs = get(gca, 'yticklabels');
%     yticklabs(3:2:end) = '';
%     set(gca, 'yticklabels', yticklabs);
% end

hax = gca;
pp = get(gcf, 'position');
if(scale_fig)
    set(gcf, 'position', [pp(1:2) 570  450]);
else
    set(gcf, 'position', [pp(1:2) 570  590]);
end



switch signal_type
  case 'bat'
    set(gcf, 'position', [pp(1:2) 570  320]);
end


for n = 1:length(hax.YTickLabels)
    if(strcmp(hax.YTickLabels{n}(1:4), 'comp'))
        % hax.YTickLabels{n} = '';
        hax.YTickLabels{n} = strrep(hax.YTickLabels{n}, 'comp. ', '');
        if((length(hax.YTickLabels) > 20) & (n > 10) & (rem(n, 2) == 0))
        % if((n ~= length(hax.YTickLabels)) & (n > 4))
            hax.YTickLabels{n} = '';
        end
    end

end
hax.XAxis.Visible = 'off';
% hax.YAxis.Visible = 'off';

set_gca_fonts(FONT_NAME, FONT_SIZE, hax);

% set colours:
llcube = cubehelix(10);
llcube = llcube(1:8, :);
% lcube = cubehelix(8, 2.5, -1.5, 3.6, 0.7, [0.2, 1], [0.2, 0.6]);
% lcube = cubehelix(16, 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]);
% lcube = (lcube(2:end - 1, :));
for n = 2:length(hp)
    in = mod(n - 2, size(llcube, 1)) + 1;
    hp(n).Color = llcube(in, :);
end

lc = cubehelix(8, 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]);
lblue = lc(2, :);
lred = lc(3, :);
hp(1).Color = lblue;

% ''components'' text
ys = ylim;
yset = (ys(2) - ys(1)) / 2;
% h = text(-42.78, mean(ys) - yset * 0.2, 'components', 'horizontalalignment', 'right', ...
%          'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
h = text(-25, mean(ys) - yset * 0.00, 'components', 'horizontalalignment', 'right', ...
         'fontname', FONT_NAME, 'fontsize', FONT_SIZE);

set(h, 'rotation', 90);


if(print_)
    fname = ['pics/' signal_type '_test/' signal_type '_comps_' method_str];
    % print([fname '_v2.svg'], '-dsvg');
    print2eps([fname '_v2.eps']);
end


%---------------------------------------------------------------------
% error plot
%---------------------------------------------------------------------

subplot = @(n, m, p) subtightplot(n, m, p, [0.1, 0.1], [0.21, 0.05], [0.1, 0.01]);

set_figure(901);
pp = get(gcf, 'position');
set(gcf, 'position', [pp(1:2) 570 296]);

hx(1) = subplot(2, 1, 1); hold all;
plot(x, '-o', 'markersize', 2, 'color', lblue);
plot(y, '-+', 'markersize', 2, 'color', lred);
ys = ylim();
hg = legend({'signal', 'estimate'}, 'location', 'southeast', ...
            'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
hg.Position = [0.76 0.555 0.2140 0.1436];

xlim([0 length(x)]);
set(hx(1), 'xticklabels', []);
hx(1).XAxis.Visible = 'off';
hx(1).Clipping = 'off';

switch signal_type
  case 'bat'
    ylim([-0.2 0.2]);
  case 'nnlfm4'
    ylim([-3.5, 3.5]);
end
ys = ylim;

hx(2) = subplot(2, 1, 2); hold all;
switch signal_type
  case 'nnlfm4'
    d = load('data/test_signals/ffgn_1_02_1_512_0_signal.mat');
    n = d.x(1:256);
    n = n .* (0.6081);


    hn = plot(n, '-+', 'markersize', 2, 'color', lred);    
    he = plot(x - y, '-o', 'markersize', 2, 'color', lblue);    


    hleg = legend([he, hn], {'residual', 'noise'}, 'location', 'southeast', ...
                  'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
    hleg.Position = [0.78 0.1769 0.1677 0.1368];
    set(hx, 'ytick', [-3.5 0 3.5]);
    
  otherwise
    plot(x - y, '-', 'color', lblue);

    legend({'residual'}, 'location', 'southeast', ...
           'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
end
hx(2).Clipping = 'off';

ylim(ys);
xlim([0 length(x)]);

set_gca_fonts(FONT_NAME, FONT_SIZE - 2, hx);

xlabel('samples', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);


if(print_)
    fname = ['pics/' signal_type '_test/' signal_type '_terror_' method_str];
    % print([fname '_v2.svg'], '-dsvg');
    print2eps([fname '_v2.eps']);
end


%---------------------------------------------------------------------
% plot TFDs for each component?
%---------------------------------------------------------------------
switch signal_type
  case 'bat'

    ggap = 0.03;
    subplot = @(n, m, p) subtightplot(n, m, p, [0.02, 0.02], [ggap, ggap], [ggap + 0.02, ggap + 0.02]);

    N = length(x);
    L = length(y_comps);
    if(L > 4)
        L = 4;
    end
    
    set_figure(903);
    pp = get(gcf, 'position');
    set(gcf, 'position', [pp(1:2) 570 200]);
    cc = cubehelix(256);
    cc = flipud(cc);
    colormap(cc); % , 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]));
    for n = 1:L
        qtfd = qtfd_sep_kern(y_comps{n}, {63, 'hann'}, {63, 'dolph', 100}, N, N);
        hx = subplot(1, L, n);
        imagesc(qtfd);
        axis('xy'); axis('tight'); axis('square');
        % hx.XAxis.Visible = 'off';
        % hx.YAxis.Visible = 'off';
        hx.Box = 'on';
        hx.XAxis.TickValues = [];
        hx.YAxis.TickValues = [];
        if(n == 1)
            xlabel('frequency', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
            ylabel('time', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
        end
        if(n == 1)
            title_right_align(['components: ' num2str(n)], FONT_NAME, FONT_SIZE);
        else
            title_right_align(num2str(n), FONT_NAME, FONT_SIZE);
        end
    end


    if(print_)
        fname = ['pics/' signal_type '_test/' signal_type '_TFDcomponents_' method_str];
        print([fname '_v2.svg'], '-dsvg');
        print2eps([fname '_v2.eps']);
    end
    

end

end

