%-------------------------------------------------------------------------------
% plot_2tone_tests_compare: load results from .mat files and plot in 1 figure
%
% Syntax: [] = plot_2tone_tests_compare()
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
% Started: 17-09-2021
%
% last update: Time-stamp: <2023-07-06 07:30:26 (otoolej)>
%-------------------------------------------------------------------------------
function plot_2tone_tests_compare(time_mask_tukey, print_)
if(nargin < 1 || isempty(time_mask_tukey)), time_mask_tukey = true; end
if(nargin < 2 || isempty(print_)), print_ = false; end



methods = {'vmd', 'tvemd', 'ssst', 'xtfd', 'tvfilt'};
methods = {'tvemd', 'vmd', 'efd', 'ssst', 'xtfd', 'tvfilt'};
    
if(time_mask_tukey)
    tstr = '_septime_tukey';            
end
else
    tstr = '';
end

FONT_NAME = 'helvetica';
FONT_SIZE = 14;


set_figure(40);
pp = get(gcf, 'position');
% set(gcf, 'position', [pp(1:2) 1000 560]);
set(gcf, 'position', [pp(1:2) 600 800]);

cc = cubehelix(256);
cc = [0.9312692223325372, 0.8201921796082118, 0.7971480974663592;
 0.909228119877378, 0.7604574726010956, 0.7498481961062412; 
 0.8845085525304266, 0.7007680763318311, 0.7096101836718632; 
 0.8559578605899612, 0.6418993116910497, 0.6754191211563135; 
 0.8226766206933027, 0.584539328313147, 0.646036981071441; 
 0.7840440880599453, 0.5292660544265891, 0.6200568926941761; 
 0.739734329496642, 0.4765280683170713, 0.5959617419736206; 
 0.689722696893931, 0.4266298828270737, 0.5721854455518277; 
 0.6390980743608128, 0.3835159423260408, 0.5493432793563523; 
 0.57916573903086, 0.33934576125314425, 0.5219003947563425; 
 0.5151069036855755, 0.29801047535056074, 0.49050619139300705; 
 0.4479571835317905, 0.25920242379465147, 0.4539864481633419; 
 0.37894937987024996, 0.2224702044652721, 0.41140014301575434; 
 0.3094687739404407, 0.18723350849627401, 0.3620818115057747; 
 0.24100281569125206, 0.15280235523757413, 0.30567558960763097; 
      0.1750865648952205, 0.11840023306916837, 0.24215989137836502];

% cc = [0.7788013041060939, 0.8977179643101579, 0.7971480974663592;
%  0.7067610774608281, 0.8634064504727101, 0.7498481961062412;
%  0.6363354543499341, 0.8269573389790794, 0.7096101836718632;
%  0.5688299932818057, 0.7878960136975732, 0.6754191211563135;
%  0.5053306243197957, 0.7459011261068167, 0.646036981071441;
%  0.44667003488651325, 0.7008115737370916, 0.6200568926941761;
%  0.39340270878404693, 0.6526282668058052, 0.5959617419736206;
%  0.34578922187214445, 0.6015106895337712, 0.5721854455518277;
%  0.30708191457572354, 0.5523371187445275, 0.5493432793563523;
%  0.2699414955609766, 0.4965778688514405, 0.5219003947563425;
%  0.23752283590185305, 0.4391544151774967, 0.49050619139300705;
%  0.2089672147705885, 0.3807223151896951, 0.4539864481633419;
%  0.18318359421158437, 0.32201177575137396, 0.41140014301575434;
%  0.15888694405339765, 0.26380026863600353, 0.3620818115057747;
%  0.13464342464737084, 0.20688320940685445, 0.30567558960763097;
%  0.10892122810595439, 0.15204350513797243, 0.24215989137836502];


cc = cubehelix(256);
% cc = cubehelix(256, 2.5, -1.5, 3.6, 0.7, [0.2, 1], [0.2, 0.6]);
% cc = flipud(cc);
colormap(cc); % , 1.5, 3, 4, 1, [0.2, 1], [0, 0.9]));


subplot = @(n, m, p) subtightplot(n, m, p, [0.04, 0.02], [0.1, 0.04], [0.12, 0.15]); % , marg_h,marg_w);


for n = 1:length(methods)
        % hxs = subplot(2, 3, n); hold all;
        hxs = subplot(3, 2, n); hold all;

        % if(strcmp(methods{n}, 'xtfd'))
        %     d = load(['./data/tones_test_' methods{n} tstr '_v3.mat']);
        %     dispVars(methods{n}, 'here');
        % else
        d = load(['./data/plots/tones_test_' methods{n} tstr '.mat']);
            % end
            fprintf(col_str('loading parameters from file saved on %s\n', 1), ...
                    datestr(d.time_now));

            rex = d.re;
            rex(rex > 0.5) = 0.5;
            imagesc(log10(d.amp_ratio), d.f_ratio, rex, [0 0.5]);
            axis('xy');
            axis('tight');
            set(hxs, 'TickDir', 'out');
            set(hxs, 'ytick', 0.2:0.2:0.8);
            set(hxs, 'fontname', FONT_NAME, 'fontsize', FONT_SIZE - 2);
            if(n > 4)
        xlabel('amplitude ratio (log)', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
        set(hxs, 'xtick', [-2:2]);
    else
        set(hxs, 'xtick', [-2:2]);
        set(hxs, 'xticklabels', []);
    end
    if(n == 1 | n == 3 | n == 5)
        ylabel('frequency ratio', 'fontname', FONT_NAME, 'fontsize', FONT_SIZE);
    else
        set(hxs, 'yticklabels', []);
    end
    xlim([-2, 2]);

    % put a box around the plots:
    do_box = true;
    if(do_box)
        xs = xlim;
        ys = ylim;
        hrec = rectangle('position', [xs(1), ys(1), xs(2) - xs(1), ys(2) - ys(1)], ...
                         'edgecolor', [0 0 0]);
        hrec.LineWidth = 0.2;
    end
    

    title(upper(strrep(methods{n}, 'tv', 'tv-')), 'fontname', FONT_NAME, 'fontsize', FONT_SIZE); % );


    % colourbar at the end
    if(n == length(methods))
        pp = get(hxs, 'position');

        hcol = colorbar('fontname', FONT_NAME, 'fontsize', FONT_SIZE - 2);
        % hcol.Layout.Tile = 'east';
        hcol.Label.String = 'relative error';
        hcol.FontSize = FONT_SIZE - 1;

        set(hxs, 'position', pp);

        % hcol.Position = [0.9288 0.1000 0.0198 0.3893];
        % hcol.Position = [0.915 0.1000 0.0198 0.3893];        
        hcol.Position = [0.865 0.1000 0.029 0.26];                
    end
end

if(print_)
    % print(['pics/two_tone_test_6methods_' tstr '_v2.svg'], '-dsvg');
    print2eps(['pics/two_tone_test_6methods_' tstr '_v5.eps']);
end


end

