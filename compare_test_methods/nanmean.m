function m = nanmean(x)
m = sum(x(~isnan(x))) ./ length(x);
