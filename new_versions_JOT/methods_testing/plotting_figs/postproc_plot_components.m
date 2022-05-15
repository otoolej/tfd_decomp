%-------------------------------------------------------------------------------
% postproc_plot_components: 
%
% Syntax:  = postproc_plot_components(x, y_comps)
%
% Inputs: 
%     x, y_comps - 
%
% Outputs: 
%      - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 11-02-2022
%
% last update: Time-stamp: <2022-04-24 14:10:15 (otoolej)>
%-------------------------------------------------------------------------------
function postproc_plot_components(x, y_comps, fignum, scale_fig)
if(nargin < 3 || isempty(fignum)), fignum = 800; end
if(nargin < 4 || isempty(scale_fig)), scale_fig = false; end


FONT_NAME = 'helvetica';
FONT_SIZE = 14;
if(scale_fig)
    FONT_SIZE = 10;
end


hp = plot_components(x, y_comps, 1, fignum, length(y_comps), false);

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
% just add a 'C' in front of components:
for p = 1:length(hax.YAxis.TickLabels) - 1
    hax.YAxis.TickLabels{p + 1} = ['C' hax.YAxis.TickLabels{p + 1}];
end


hax.XAxis.Visible = 'off';
set_gca_fonts(FONT_NAME, FONT_SIZE, hax);

% set colours:
lcube = cubehelix(10);
lcube = lcube(1:8, :);
% lcube = cubehelix(8, 2.5, -1.5, 3.6, 0.7, [0.2, 1], [0.2, 0.6]);
% lcube = cubehelix(16, 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]);
% lcube = (lcube(2:end - 1, :));
for n = 2:length(hp)
    in = mod((2 * n) - 2, size(lcube, 1)) + 1;
    hp(n).Color = lcube(in, :);
end

lc = cubehelix(8, 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]);
lblue = lc(2, :);
lred = lc(3, :);
hp(1).Color = lblue;

% ''components'' text
do_comp_txt = false;
if(do_comp_txt)
    ys = ylim;
    yset = (ys(2) - ys(1)) / 2;
    % h = text(-42.78, mean(ys) - yset * 0.2, 'components', 'horizontalalignment', 'right', ...
    %          'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
    h = text(-105, mean(ys) - yset * 0.2, 'components', 'horizontalalignment', 'right', ...
             'fontname', FONT_NAME, 'fontsize', FONT_SIZE);

    set(h, 'rotation', 90);
end




