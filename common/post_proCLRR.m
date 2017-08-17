function grp = post_proCLRR(C,K,d,alpha)
warning off
if(nargin<4)
	alpha = 4;
end
if(nargin<3)
	d = 4;
end
%post processing
[U,S,V] = svd(C,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^4;
% spectral clustering
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,S,V] = svd(L);
V = U(:,1:K);
V = D*V;
grp = kmeans(V,K,'emptyaction','singleton','replicates',20,'display','off');