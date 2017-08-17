% roml using feature only
% F is a 1xF cell containing dxN sift features
% N is the number of inlier points
function [PM,D,E,L] = roml_feature(F,N,lambda,P_init)
d = size(F{1},1);
K = size(F,2);
if(nargin<2)
	N = size(F{1},2);
end
if(nargin<3)
	lambda = 5/sqrt(N);
end
if(nargin<4)
	for i=1:K
		P_init{i} = speye(N);
	end
end


X = []; %inlier feature matrix: dNxK
for i = 1:K
	Nk = size(F{i},2);
	PM{i} = speye(Nk,N);
	X = [X vec(F{i}*P_init{i})];
	G{i} = sparse(kron(speye(N),F{i})); % Kronecker product
	J{i} = sparse(kron(speye(N),ones(1,Nk)));
	H{i} = sparse(kron(ones(1,N),speye(Nk)));
end

% initiate parameters for ADMM
%rho = 1e-8; 
rho = 1e-2;
eta = 1.01;
max_rho = 1e10; maxIter = 1e8;
epsilon = 1e-6;

IK = eye(K);
D = X;
L = zeros(d*N,K);
E = zeros(d*N,K);
Y = zeros(d*N,K);

% start main loop
iter = 0;
disp(['initial,rank = ' num2str(rank(L))]);
while(iter<maxIter)
	iter = iter + 1;
	% update L
	temp = D-E+Y/rho;
	[U,sigma,V] = svd(temp,'econ');
	sigma = diag(sigma);
	svp = length(find(sigma>1/rho));
	if(svp>=1)
		sigma = sigma(1:svp)-1/rho;
	else
		svp = 1;
		sigma = 0;
	end
	L = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
	% update E
	temp = D-L+Y/rho;
	E = max(0,temp-lambda/rho)+min(0,temp+lambda/rho);
	%E = sign(temp).*max(abs(temp)-lambda/rho,0);
	% update Pk
	temp = [];
	%P1 = speye(size(F{1},2),N);
	%PM{1} = P1;
	%temp = [temp G{1}*P1(:)];
		
	for k = 1:K
		Nk = size(F{k},2);		
		f = IK(k,:)*(Y'-rho*(L+E)')*G{k};
		
		fMat = reshape(f,Nk,N);
		fMat = [fMat,Inf*ones(Nk,N-Nk)];
		cost = fMat+abs(min(fMat(:)));		
		assignment = assignmentoptimal(cost);		
		assignment = assignment(1:N,:);
		INk = speye(Nk);
		theta = INk(assignment,:);        
		temp = [temp G{k}*theta(:)];
		PM{k} = theta;
		%toc
	end	
		
	D = temp;
	
	leq = D-L-E;
	stopC = max(max(abs(leq)));
	if(iter==1 || mod(iter,50)==0 || stopC<epsilon)
		disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ...
			',rank=' num2str(rank(L,1e-3*norm(L,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
	end
	if(stopC<epsilon)
		break;
	else
		Y = Y + rho*leq;
		rho = min(max_rho,rho*eta);
	end
end
test = 1;