function [V,lambda]=compute_embedding(W,Rdim)

D=diag(sum(W,2));
L=D-W;
L(isnan(L)) = 0; D(isnan(D)) = 0;
L(isinf(L)) = 0; D(isinf(D)) = 0;
tol=0;
options.disp = 0;
options.isreal = 1;
options.issym = 1;
clear W;
[V, lambda] = eigs(L, D, Rdim + 1,0,options);

% Sort eigenvectors in ascending order
lambda = diag(lambda);
[lambda, ind] = sort(lambda, 'ascend');
lambda = lambda(2:Rdim + 1);

% Final embedding
V = V(:,ind(2:Rdim + 1));