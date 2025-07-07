function [y, breaks_new, breaks_lable] = Transform_for_uneven_colorbar_percentile(x, ...
    breaks_lable_round)

% y is the transformed data matrix used to draw the map, sort index

[~, I] = sort(x(:)); 
[~, P] = sort(I); % position
y = reshape(P, size(x));
y(isnan(x)) = NaN;
%ticks = 0:0.1:1;
ticks = [0:0.1:0.9, 0.9975, 1];
breaks_new = quantile(y(:), ticks);
breaks_lable = quantile(x(:), ticks);
breaks_lable = string(round(breaks_lable, breaks_lable_round, "significant"));
breaks_lable(1) = " ";
breaks_lable(end) = " ";
end