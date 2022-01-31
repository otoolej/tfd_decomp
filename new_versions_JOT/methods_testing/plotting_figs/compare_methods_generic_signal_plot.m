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
% last update: Time-stamp: <2021-12-06 10:06:01 (otoolej)>
%-------------------------------------------------------------------------------
function compare_methods_generic_signal_plot(method_str, signal_type, print_)
if(nargin < 1 || isempty(method_str)), method_str = 'vmd'; end
if(nargin < 2 || isempty(signal_type)), signal_type = 'bat'; end
if(nargin < 3 || isempty(print_)), print_ = false; end



FONT_NAME = 'helvetica';
FONT_SIZE = 12;

all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'msst', 'vmd'};

% [x, x_components, y, y_comps] = compare_methods_testsignals('noise', {'tvfilt'});
[x, x_components, y, y_comps] = compare_methods_testsignals(signal_type, {method_str});

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
set(gcf, 'position', [pp(1:2) 570  590]);

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
lcube = cubehelix(10);
lcube = lcube(1:8, :);
% lcube = cubehelix(8, 2.5, -1.5, 3.6, 0.7, [0.2, 1], [0.2, 0.6]);
% lcube = cubehelix(16, 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]);
% lcube = (lcube(2:end - 1, :));
for n = 2:length(hp)
    in = mod(n - 2, size(lcube, 1)) + 1;
    hp(n).Color = lcube(in, :);
end

lc = cubehelix(8, 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]);
lblue = lc(2, :);
lred = lc(3, :);
hp(1).Color = lblue;

% ''components'' text
ys = ylim;
yset = (ys(2) - ys(1)) / 2;
h = text(-42.78, mean(ys) - yset * 0.2, 'components', 'horizontalalignment', 'right', ...
         'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
set(h, 'rotation', 90);


if(print_)
    fname = ['pics/' signal_type '_test/' signal_type '_comps_' method_str];
    print([fname '.svg'], '-dsvg');
    print2eps([fname '.eps']);
end


%---------------------------------------------------------------------
% error plot
%---------------------------------------------------------------------

subplot = @(n, m, p) subtightplot(n, m, p, [0.1, 0.1], [0.2, 0.05], [0.1, 0.01]);

set_figure(901);
pp = get(gcf, 'position');
set(gcf, 'position', [pp(1:2) 570 296]);

hx(1) = subplot(2, 1, 1); hold all;
plot(x, '-o', 'markersize', 2, 'color', lblue);
plot(y, '-+', 'markersize', 2, 'color', lred);
ys = ylim();
hg = legend({'signal', 'estimate'}, 'location', 'southeast', ...
            'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
hg.Position = [0.7760 0.5295 0.2140 0.1436];

xlim([0 length(x)]);
set(hx(1), 'xticklabels', []);
hx(1).XAxis.Visible = 'off';

switch signal_type
  case 'bat'
    ylim([-0.2 0.2]);
end

hx(2) = subplot(2, 1, 2); hold all;
plot(x - y, '-', 'color', lblue);
ylim(ys);
legend({'error'}, 'location', 'southeast', ...
       'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
xlim([0 length(x)]);
xlabel('samples');

set_gca_fonts(FONT_NAME, FONT_SIZE, hx);




if(print_)
    fname = ['pics/' signal_type '_test/' signal_type '_terror_' method_str];
    print([fname '.svg'], '-dsvg');
    print2eps([fname '.eps']);
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
        print([fname '.svg'], '-dsvg');
        print2eps([fname '.eps']);
    end
    

end
