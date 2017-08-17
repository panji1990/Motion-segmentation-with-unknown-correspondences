% This function solves the following problem:
% min_{C,Pf,L,M} 1/2||C||_F^2 + lambda1||E1||_1 +
% lambda2||M||_*+lambda3||E2||_1 s.t., D1=L+E1, L=LC, D2=M+E2
function [C,P,L,M,D1,D2,E1,E2] = JEPS_feature(XC,XF,lambda1,lambda2,lambda3,affine,P_init)
[d1,N] = size(XC);
d2 = size(XF{1},1);
F = d1/2;
if(nargin<3)
	lambda1 = 5/sqrt(N);
end
if(nargin<4)
	lambda2 = 1;
end
if(nargin<5)
	lambda3 = lambda1;
end
if(nargin<6)
	affine = true;
end
if(nargin<7)
	for i=1:F
		P_init{i} = speye(N);
	end
end

% initial parameters for ADMM
rho = 1e-2; 
eta = 1.01;
max_rho = 1e10; maxIter = 1e8;
epsilon = 1e-6;

IN = eye(N);
Id1 = eye(d1);
IF = eye(F);
ONES = ones(N,N);

D1 = XC;
D2 = [];
for i=1:F
	D2 = [D2 reshape(XF{i}*P_init{i},d2*N,1)];
	D1((2*i-1):(2*i),:) = XC((2*i-1):(2*i),:)*P_init{i};
	G{i} = sparse(kron(speye(N),XF{i}));
end
J = sparse(kron(speye(N),ones(1,N)));
H = sparse(kron(ones(1,N),speye(N)));

L = zeros(d1,N);
M = zeros(d2*N,F);
C = zeros(N,N);
E1 = zeros(d1,N);
E2 = zeros(d2*N,F);
Y1 = zeros(d1,N);
Y2 = zeros(d1,N);
Y3 = zeros(d2*N,F);
Y4 = zeros(1,N);

% start main loop
iter = 0;
disp(['initial,rank of L = ' num2str(rank(L)) ', rank of M = ' num2str(rank(M))]);
while(iter<maxIter)
	iter=iter+1;	
	%update L
	%lhs = lambda1*IN+rho*(IN-C)*(IN-C');
	%rhs = lambda1*D1-Y1+Y1*C';
	lhs = rho*(IN+(IN-C)*(IN-C'));
	rhs = Y1-Y2*(IN-C')+rho*(D1-E1);
	L = rhs/lhs;	
	%update C
	rltl = rho*(L'*L+affine*ONES);
	lhs = IN+rltl;
	rhs = L'*Y2+rltl-affine*ones(N,1)*Y4;
	C = lhs\rhs;
% 	[~,sigma,V] = svd(L,'econ');
% 	sigma = diag(sigma);
% 	svp = max(1,length(find(sigma>1e-1*sigma(1))));	
% 	C = V(:,1:svp)*V(:,1:svp)';	
	%update M
	temp = D2-E2+Y3/rho;
	try
	   [U,sigma,V] = svd(temp,'econ');
	catch err
	   [U,sigma,V] = svd((temp+temp')/2,'econ');
	end
	sigma = diag(sigma);
	svp = length(find(sigma>lambda2/rho));
	if(svp>=1)
		sigma = sigma(1:svp)-lambda2/rho;
	else
		svp = 1;
		sigma = 0;
	end
	M = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
	%update E1
	temp = D1-L+Y1/rho;
	E1 = max(0,temp-lambda1/rho)+min(0,temp+lambda1/rho);
	%update E2
	temp = D2-M+Y3/rho;
	E2 = max(0,temp-lambda3/rho)+min(0,temp+lambda3/rho);
	%update Pf
	D2 = [];
	%D2 = vec(XF{1});
	%P{1} = sparse(IN);
	for k=1:F
		alfa = Id1(:,2*k-1:2*k);
		W = XC(2*k-1:2*k,:);
		%cost1 = -lambda1*(L'*alfa*W)';
		cost1 = ((Y1'-rho*(L'+E1'))*alfa*W)';
		e = IF(k,:);
		% 		%cost2 = e*(Y2'-rho*(M+E)')*G{k};
		cost2 = e*(Y3'-rho*(M+E2)')*G{k};
		cost2 = reshape(cost2,N,N);
		cost = cost1+cost2;
		cost_pos = cost+abs(min(cost(:)));
		assignment = assignmentoptimal(cost_pos);
		%assignment = assignmentsuboptimal1(cost_pos);
 		% f = cost1(:)'+cost2;
        % option = optimset('LargeScale','on','Display','off');
 		% theta = linprog(f,H,ones(N,1),J,ones(N,1),zeros(N*N,1),[],[],option);
 		% theta = round(theta);
 		% perm = sparse(reshape(theta,N,N));
		perm = sparse(IN(assignment,:));
		P{k} = perm;
		D1((2*k-1):(2*k),:) = XC((2*k-1):(2*k),:)*perm; %update D1
		D2 = [D2 reshape(XF{k}*perm,d2*N,1)]; %update D2
	end
	
	leq1 = D1-L-E1;
	leq2 = L-L*C;
	leq3 = D2-M-E2;	
	leq4 = affine*(ones(1,N)*C-ones(1,N));
	stopC1 = max(norm(abs(leq1(:)),Inf),norm(abs(leq2(:)),Inf));
	stopC2 = max(norm(abs(leq3(:)),Inf),norm(abs(leq4(:)),Inf));
	stopC = max(stopC1,stopC2);
	if(iter==1 || mod(iter,50)==0 || stopC<epsilon)
		disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ...
			',rankL=' num2str(rank(L,1e-3*norm(L,2))) ',rankM='...
			num2str(rank(M,1e-3*norm(M,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
	end
	if(stopC<epsilon)
		break;
	else
		Y1 = Y1+rho*leq1;
		Y2 = Y2+rho*leq2;
		Y3 = Y3+rho*leq3;
		Y4 = Y4+rho*leq4;
		rho = min(max_rho,rho*eta);
	end
end










