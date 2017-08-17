function Q1trim= comp_affine_invariant_coord(X)
n=length(X);
PA=X*inv((X'*X))*X';
[Q1,R1, E1] = qr(PA);
Q1trim=Q1(:,1:3);

