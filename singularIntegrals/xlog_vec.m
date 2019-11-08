function [ out ] = xlog_vec(v)
% returns the vector v log(|v|), and 0 if |v| = 0
% if v is a matrix, the norm is taken across the second dimension so that
% the operation is applied column-wise. 

norm_v = repmat(sqrt((v(:,1).^2 + v(:,2).^2)),1,2);
out = v.*log(norm_v);
out(isnan(out)) = 0;




end

