function [MissrateNC,C,grpNC] = lrr_El2_cvx(X,s,lambda,opt)
if(nargin<4)
	opt = 'Noisy';
end
if(nargin<3)
	lambda = 0.01;
end

[d,N] = size(X);

if(strcmp(opt,'Noisy'))
	cvx_begin;
	cvx_precision high
	variable C(N,N)
	minimize( norm_nuc(C) + lambda*norm(X-X*C,'fro') );
	subject to
	sum(C) == 1;
	%diag(C) == 0;
	C(:) >= 0;
	cvx_end;
elseif(strcmp(opt,'Noiseless'))
	cvx_begin;
	cvx_precision high
	variable C(N,N)
	minimize( trace(C)  );
	subject to
	X == X*C;
	%sum(C) == 1;
	%diag(C) == 0;
	C(:) >= 0;
	cvx_end;
else
	disp('Two options of opt: Noisy or Noiseless');
end
C = full(C);
W = abs(C) + abs(C');
K = max(s);
grpNC = SpectralClustering(W,K);
MissrateNC = ErrorRate(grpRC,s);