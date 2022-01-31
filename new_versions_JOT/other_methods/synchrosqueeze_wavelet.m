%-------------------------------------------------------------------------------
% synchrosqueeze_wavelet: decomposition using the inverse wavelet synchrosqueeze
% transform
%
% Syntax: [y, x_res, x_components] = synchrosqueeze_wavelet(x, Fs, N_components, db_plot)
%
% Inputs: 
%     x, Fs, N_components, db_plot - 
%
% Outputs: 
%     [y, x_res, x_components] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 02-09-2021
%
% last update: Time-stamp: <2021-09-09 17:17:10 (otoolej)>
%-------------------------------------------------------------------------------
function [y, x_res, y_comps] = synchrosqueeze_wavelet(x, Fs, N_components, db_plot, wavelet_type, ridge_pen)
if(nargin < 2 || isempty(Fs)), Fs = 1; end
if(nargin < 3 || isempty(N_components)), N_components = 1; end
if(nargin < 4 || isempty(db_plot)), db_plot = false; end
if(nargin < 5 || isempty(wavelet_type)), wavelet_type = 'amor'; end
if(nargin < 6 || isempty(ridge_pen)), ridge_pen = 0; end



x = x(:)';



% wavelet synchrosqueeze transform:
[tf_swav, freq] = wsst(x, Fs, wavelet_type);

% estimate ridges in TF plane:
[fridge, iridge] = wsstridge(tf_swav, freq, ridge_pen, 'NumRidges', N_components);

% invert 
y_comps = iwsst(tf_swav, iridge, wavelet_type);

% add components
y = nansum(y_comps, 2)';
x_res = x - y;



y_comps = num2cell(y_comps, 1);

if(db_plot)
    set_figure(202);
    t = (1:length(x)) ./ Fs;
    pcolor(t, freq, abs(tf_swav));
    shading interp;
    axis('tight');
    % grid on;
    title('Wavelet Synchrosqueezed Transform');
    ylabel('Frequency');
    xlabel('Time');

    plot(t, fridge,'k','linewidth',2);
    view(90, -90);

    set_figure(201);
    subplot(2, 1, 1); hold all;    
    title('Wavelet Synchrosqueezed Transform');        
    plot(x, '-o');
    plot(y, '-+');
    ys = ylim();

    subplot(2, 1, 2); hold all;
    plot(x - y, '-');
    ylim(ys);

    fprintf('MSE = %g\n', nanmean(abs(x - y) .^ 2));

    plot_components(x, y_comps, Fs, 25, N_components, false);
    title('Wavelet Synchrosqueezed Transform');            
end

