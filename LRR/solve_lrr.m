% inputs: X -- data matrix
%         lambda -- weight for noise term E
%         reg_type -- regularization type of E
%                  -- reg = 0 (default),    use L21-norm
%                  -- reg = 1 (or nonzero), use L1-norm
%         affine -- affine constraint
%         r -- projection dimension
%           -- r = 0(default), no projection
%           -- r = nonzero, project with PCA
%         post_pro -- post processing C
% outputs: Missrate -- missclassification rate
%          C -- coefficient matrix
%          E -- noise matrix
%
% Pan Ji, pan.ji@anu.edu.au
% Sep 2013, ANU

function [Missrate,C,E,grp,time] = solve_lrr(X,s,lambda,reg_type,affine,r,post_pro)
if(nargin < 7 || isempty(post_pro))
    post_pro = 0;
end

if(nargin < 6 || isempty(r))
    r = 0;
end

if(nargin < 5 || isempty(affine))
    affine = 0;
end

if(nargin < 4 || isempty(reg_type))
    reg_type = 0;
end

if(nargin < 3 || isempty(lambda))
    lambda = 0.24;
end

Xp = DataProjection(X,r);

% Q = orth(Xp');
% A = Xp*Q;

if(reg_type == 0)
    if(affine)
		tic
        [C,E] = alm_lrr_l21_affine(Xp,Xp,lambda);
		time = toc;
	else
		tic
        [C,E] = alm_lrr_l21(Xp,Xp,lambda);
		time = toc;
    end
elseif(reg_type == 1)
    if(affine)
        [C,E] = alm_lrr_l1_affine(Xp,Xp,lambda);
    else
        [C,E] = alm_lrr_l1(Xp,Xp,lambda);
	end
elseif(reg_type == 2)
	[C,E] = alm_lrr_l2_affine(Xp,Xp,lambda);
end

% C = Q*C;
K = max(s);
if(post_pro == 0)
    W = abs(C)+abs(C');	
	grp = SpectralClustering(W,max(s));
	
else    
    % refining C   
	try
		[U,S,V] = svd(C,'econ');
	catch err
		[U,S,V] = svd((C+C')/2,'econ');
	end
    S = diag(S);
    r = min(4*K+1,sum(S>1e-3*S(1)));
    S = S(1:r);
    U = U(:,1:r)*diag(sqrt(S));
    U = normr(U);
    Z = U*U';Z=abs(Z);L = Z.^4;
	    
    % spectral clustering
    L = (L + L')/2;
    D = diag(1./sqrt(sum(L,2)));
    L = D*L*D;
    [U,S,V] = svd(L,'econ');
    
    V = U(:,1:K);
    V = D*V;
    grp = kmeans(V,K,'emptyaction','singleton','replicates',20,'display','off');
end

%--------------------------------------------------------------------------
Missrate = ErrorRate(grp,s);
end