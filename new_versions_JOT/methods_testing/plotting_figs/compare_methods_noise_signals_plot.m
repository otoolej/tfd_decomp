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
% last update: Time-stamp: <2021-11-19 12:35:40 (otoolej)>
%-------------------------------------------------------------------------------
function compare_methods_noise_signals_plot(method_str, print_)
if(nargin < 1 || isempty(method_str)), method_str = 'vmd'; end
if(nargin < 2 || isempty(print_)), print_ = false; end



FONT_NAME = 'helvetica';
FONT_SIZE = 10;

all_methods = {'tvfilt', 'xtfd', 'efd', 'tvemd', 'ssst', 'msst', 'vmd'};

% [x, x_components, y, y_comps] = compare_methods_testsignals('noise', {'tvfilt'});
[x, x_components, y, y_comps] = compare_methods_testsignals('noise', {method_str});

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

for n = 1:length(hax.YTickLabels)
    if(strcmp(hax.YTickLabels{n}(1:4), 'comp'))
        % hax.YTickLabels{n} = '';
        hax.YTickLabels{n} = strrep(hax.YTickLabels{n}, 'comp. ', '');
        if((length(hax.YTickLabels) > 20) & (n > 10) & (rem(n, 2) == 0))
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

if(print_)
    print(['pics/noise_test/noise_ffgn_test_comps_' method_str '.svg'], '-dsvg');
    print2eps(['pics/noise_test/noise_ffgn_test_comps_' method_str '.eps']);
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
legend({'signal', 'estimate'}, 'location', 'southeast', ...
      'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
xlim([0 length(x)]);
set(hx(1), 'xticklabels', []);
hx(1).XAxis.Visible = 'off';

hx(2) = subplot(2, 1, 2); hold all;
plot(x - y, '-', 'color', lblue);
ylim(ys);
legend({'error'}, 'location', 'southeast', ...
       'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
xlim([0 length(x)]);
xlabel('samples');

set_gca_fonts(FONT_NAME, FONT_SIZE, hx);




if(print_)
    print(['pics/noise_test/noise_ffgn_test_error_' method_str '.svg'], '-dsvg');
    print2eps(['pics/noise_test/noise_ffgn_test_error_' method_str '.eps']);
end
