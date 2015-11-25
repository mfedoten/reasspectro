function X = reassignment_crop_matrix(X,percent)

lim = prctile(X(:),percent);
max_val = max(X(X(:)<lim));
X(X>max_val) = max_val;