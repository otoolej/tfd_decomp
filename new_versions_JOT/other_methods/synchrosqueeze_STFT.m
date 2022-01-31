%-------------------------------------------------------------------------------
% synchrosqueeze_STFT: decomposition using the inverst STFT synchrosqueeze transform
%
% Syntax: [y, x_res, y_comps] = synchrosqueeze_STFT(x, Fs, N_components, db_plot, stft_win)
%
% Inputs: 
%     x, Fs, N_components, db_plot, stft_win - 
%
% Outputs: 
%     [y, x_res, y_comps] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 03-09-2021
%
% last update: Time-stamp: <2021-11-23 18:11:07 (otoolej)>
%-------------------------------------------------------------------------------
function [y, x_res, y_comps] = synchrosqueeze_STFT(x, Fs, N_components, db_plot, stft_win, ridge_pen, nbin_remove)
if(nargin < 2 || isempty(Fs)), Fs = 1; end
if(nargin < 3 || isempty(N_components)), N_components = 1; end
if(nargin < 4 || isempty(db_plot)), db_plot = false; end
if(nargin < 5 || isempty(stft_win)), stft_win = []; end
if(nargin < 6 || isempty(ridge_pen)), ridge_pen = 0; end
if(nargin < 7 || isempty(nbin_remove)), nbin_remove = []; end


if(N_components == 1)
    nbin_remove = [];
end



x = x(:)';



% STFT synchrosqueeze transform:
[tf_swav, freq] = fsst(x, Fs, stft_win);


% estimate ridges in TF plane:
[fridge, iridge] = tfridge(tf_swav, freq, ridge_pen, 'NumRidges', N_components, ...
                           'NumFrequencyBins', nbin_remove);

% invert 
y_comps = ifsst(tf_swav, stft_win, iridge);

% add components
y = nansum(y_comps, 2)';
x_res = x - y;



y_comps = num2cell(y_comps, 1);

if(db_plot)
    set_figure(302);
    t = (1:length(x)) ./ Fs;
    pcolor(t, freq, abs(tf_swav));
    shading interp;
    axis('tight');
    % grid on;
    title('STFT Synchrosqueezed Transform');
    ylabel('Frequency');
    xlabel('Time');

    plot(t, fridge,'k','linewidth',2);
    view(90, -90);

    set_figure(301);
    subplot(2, 1, 1); hold all;
    title('STFT Synchrosqueezed Transform');    
    plot(x, '-o');
    plot(y, '-+');
    ys = ylim();

    subplot(2, 1, 2); hold all;
    plot(x - y, '-');
    ylim(ys);

    fprintf('MSE = %g\n', nanmean(abs(x - y) .^ 2));

    plot_components(x, y_comps, Fs, 303, N_components, false);
    title('STFT Synchrosqueezed Transform');
end

