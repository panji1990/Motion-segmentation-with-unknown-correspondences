
function Z = norm_nuc( X )

%NORM_NUC   Compute nuclear norm of a matrix.

Z = sum(svd(X));
end