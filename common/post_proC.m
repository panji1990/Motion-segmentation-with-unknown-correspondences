function [grp, time] = post_proC(C,K,d,alpha)
warning off
if(nargin<4)
	alpha = 4;
end
if(nargin<3)
	d = 4;
end
tic
C = (C + C')/2;
r = d*K+1;
[U,S,~] = svd(C,'econ');
S = diag(S);
S = S(1:r);
U = U(:,1:r)*diag(sqrt(S));
U = normr(U);
Z = U*U';
L = Z.^alpha;

% ncut
L = (L + L')/2;
time = toc;
grp = ncutW(L,K); % get your ncutW funtion via http://www.cis.upenn.edu/~jshi/software/
