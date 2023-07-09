function err_plot(x, y_comps, y, fignum, method_str)
%---------------------------------------------------------------------
% plot components
%---------------------------------------------------------------------
if(nargin < 4 || isempty(fignum)), fignum = 500; end

x = x(:)';
y = y(:)';

set_figure(fignum);
subplot(2, 1, 1); hold all;
title(method_str);
plot(x, '-o');
plot(y, '-+');
ys = ylim();
xlim([0 length(x)]);

subplot(2, 1, 2); hold all;
title('error');
plot(x - y, '-');
ylim(ys);
xlim([0 length(x)]);

% fprintf('MSE = %g\n', nanmean(abs(x - y) .^ 2));

plot_components(x, y_comps, 1, fignum + 1, length(y_comps), false);
title(method_str);
