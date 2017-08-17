%--------------------------------------------------------------------------
% This is the function to call the sparse optimization program, to call the 
% spectral clustering algorithm and to compute the clustering error.
% r = projection dimension, if r = 0, then no projection
% affine = use the affine constraint if true
% s = clustering ground-truth
% missrate = clustering error
% CMat = coefficient matrix obtained by SSC
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function [missrate,CKSym,C,grps,time] = SSC(X,r,affine,alpha,outlier,rho,s)

if (nargin < 6)
    rho = 1;
end
if (nargin < 5)
    outlier = false;
end
if (nargin < 4)
    alpha = 20;
end
if (nargin < 3)
    affine = false;
end
if (nargin < 2)
    r = 0;
end

n = max(s);
Xp = DataProjection(X,r);

if (~outlier)
	tic
    CMat = admmLasso_mat_func(Xp,affine,alpha);
    C = CMat;
	time = toc;
else
	tic
    CMat = admmOutlier_mat_func(Xp,affine,alpha);
    N = size(Xp,2);
    C = CMat(1:N,:);
	time = toc;
end

%grps = 0;missrate = 1;CKSym = C;
CKSym = BuildAdjacency(thrC(C,rho));
%grps = SpectralClustering(CKSym,n);
%missrate = Misclassification(grps,s);
grps = ncutW(CKSym,n);
missrate = ErrorRate(grps,s);