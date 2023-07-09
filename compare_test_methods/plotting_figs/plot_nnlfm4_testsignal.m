%-------------------------------------------------------------------------------
% plot_nnlfm4_testsignal: 
%
% Syntax: [] = plot_nnlfm4_testsignal()
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
% Started: 15-02-2022
%
% last update: Time-stamp: <2023-07-09 02:37:16 (otoolej)>
%-------------------------------------------------------------------------------
function plot_nnlfm4_testsignal(signal_type, print_)
if(nargin < 1 || isempty(signal_type)), signal_type = 'nnlfm4'; end
if(nargin < 2 || isempty(print_)), print_ = false; end

scale_fig = true;

FONT_NAME = 'helvetica';
FONT_SIZE = 12;
if(scale_fig)
    FONT_SIZE = 15;
end



[x, y_comps, Fs, all_params] = set_signal_parameters(signal_type, false);


%---------------------------------------------------------------------
% plot:
%---------------------------------------------------------------------
y_comps = y_comps([2, 1, 3]);

switch signal_type
  case 'nnlfm4'
    x_noisy = x;
    x = set_signal_parameters('nlfm4', false);    
    d = load('data/test_signals/ffgn_1_02_1_512_0_signal.mat');
    n = d.x(1:256);
    n = n .* (0.6081);
    y_comps{4} = n;
end
    
hp = plot_components(x, y_comps, 1, 800, length(y_comps), false);

switch signal_type
  case 'nnlfm4'
    ytl = get(gca, 'yticklabels');
    ytl{5} = 'noise';
    set(gca, 'yticklabels', ytl);
    
    for n = 1:length(hp)
        hp(n).LineWidth = 1.0;
    end
end

hax = gca;
pp = get(gcf, 'position');
if(scale_fig)
    set(gcf, 'position', [pp(1:2) 570  450]);
else
    set(gcf, 'position', [pp(1:2) 570  590]);
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

set_gca_fonts(FONT_NAME, FONT_SIZE, hax);

% set colours:
llcube = cubehelix(7);
llcube = llcube(1:7, :);
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
h = text(-25, mean(ys) - yset * 0.3, 'components', 'horizontalalignment', 'right', ...
         'fontname', FONT_NAME, 'fontsize', FONT_SIZE);

set(h, 'rotation', 90);


if(print_)
    fname = ['pics/' signal_type '_test/' signal_type '_comps_justsignal_v3'];
    % print([fname '.svg'], '-dsvg');
    print2eps([fname '.eps']);
end

%---------------------------------------------------------------------
% plot TFD
%---------------------------------------------------------------------
subplot = @(n, m, p) subtightplot(n, m, p, [0.1, 0.1]);

set_figure(903);
pp = get(gcf, 'position');
set(gcf, 'position', [pp(1:2) 570 366]);
cc = cubehelix(256);
cc = flipud(cc);
colormap(cc); % , 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]));

xnn = set_signal_parameters('nlfm4', false);

hx = subplot(1, 2, 1); hold all;
N = length(x);
qtfd = qtfd_sep_kern(xnn, {63, 'hann'}, {127, 'dolph', 100}, N, N);
imagesc(qtfd);
% hx = gca;
axis('xy'); axis('tight'); axis('square');
% hx.XAxis.Visible = 'off';
% hx.YAxis.Visible = 'off';
hx.Box = 'on';
hx.XAxis.TickValues = [];
hx.YAxis.TickValues = [];

xlabel('frequency', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE + 1);
ylabel('time', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE + 1);
hx.Box = 'on';

hx = subplot(1, 2, 2); hold all;
qtfd = qtfd_sep_kern(x_noisy, {63, 'hann'}, {127, 'dolph', 100}, N, N);
imagesc(qtfd);
% hx = gca;
axis('xy'); axis('tight'); axis('square');
% hx.XAxis.Visible = 'off';
% hx.YAxis.Visible = 'off';
hx.Box = 'on';
hx.XAxis.TickValues = [];
hx.YAxis.TickValues = [];
hx.Box = 'on';


if(print_)
    fname = ['pics/' signal_type '_test/' signal_type '_TFDcomponents_justsignal_v3'];
    % print([fname '.svg'], '-dsvg');
    print2eps([fname '.eps']);
end
