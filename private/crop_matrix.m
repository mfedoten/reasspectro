function X = crop_matrix(X,percent)
% Crops large values of the input matrix according to provided percentile. The
% values above this percentile will be assigned the maximal value below the
% percentile.

lim = prctile(X(:),percent);
max_val = max(X(X(:)<lim));
X(X>max_val) = max_val;