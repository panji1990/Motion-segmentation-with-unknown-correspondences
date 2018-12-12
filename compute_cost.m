% compute cost matrix for Hungarian algorithm min \|Gvec(P)-a\|_2^2
function cost = compute_cost(G,a)
[~,N] = size(G);

tmp = (G-a*ones(1,N)).^2;
cost = sum(tmp,1); 