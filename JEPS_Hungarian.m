% This function solves the following problem:
% min_{C,Pf,L,M} 1/2||C||_F^2 + lambda1||E1||_1 +
% lambda2||M||_*+lambda3||E2||_1 s.t., D1=L+E1, L=LC, D2=M+E2
% XC, XF: Fx1 cell
function [C,P,L,M,D1,D2,E1,E2,conv] = JEPS_Hungarian(XC,XF,N,affine,param,P_init)
%[d1,N] = size(XC);
F = length(XC);
d1 = 2*F;
d2 = size(XF{1},1);
if(nargin<3)
    N = size(XF{1},2);
end
if(nargin<4)
	affine = true;
end
if(nargin<5)
	param.rho = 1e-6;
	param.eta = 1.01;
	param.epsilon = 1e-6;
	param.lambda1 = 0.05; param.lambda2 = 1; param.lambda3 = 5/sqrt(N);
end
if(nargin<6)
	for i=1:F
        Ni = size(XF{i},2);
		P_init{i} = speye(Ni,N);
	end
end

% initial parameters for ADMM
rho = param.rho; 
eta = param.eta;
epsilon = param.epsilon;
max_rho = 1e10; maxIter = 1e8;
lambda1 = param.lambda1;
lambda2 = param.lambda2;
lambda3 = param.lambda3;

IN = eye(N);
Id1 = eye(d1);
IF = eye(F);
ONES = ones(N,N);

D1 = [];
D2 = [];
for i=1:F
	D2 = [D2 reshape(XF{i}*P_init{i},d2*N,1)];
	%D1((2*i-1):(2*i),:) = XC((2*i-1):(2*i),:)*P_init{i};
    D1 = [D1; XC{i}*P_init{i}];
    G{i} = sparse(kron(speye(N),XC{i}));
	J{i} = sparse(kron(speye(N),XF{i}));
end

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
	lhs = rho*(IN+(IN-C)*(IN-C'));
	rhs = Y1-Y2*(IN-C')+rho*(D1-E1);
	L = rhs/lhs;	
	%update C
	rltl = rho*(L'*L+affine*ONES);
	lhs = IN+rltl;
	rhs = L'*Y2+rltl-affine*ones(N,1)*Y4;
	C = lhs\rhs;	
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
	D1 = [];    
%     P{1} = speye(size(XF{1},2),N);
% 	D2 = [D2 reshape(XF{1},d2*N,1)];
% 	D1 = [D1; XC{1}];
	for k=1:F
        Ni = size(XF{k},2);
        INi = [IN zeros(N,1)];
		alfa = Id1(:,2*k-1:2*k);		
		e = IF(:,k);        
        a = vec(alfa'*(L+E1-Y1/rho));        
        %b = sparse((M+E2-Y3/rho)*e);
        %tic
        cost1 = compute_cost(G{k},a);          
        %cost2 = compute_cost(J{k},b);  
        cost2 = e'*(Y3'-rho*(M'+E2'))*J{k};
        cost = reshape(cost1+cost2,Ni,N);%toc
        cost = cost-min(cost(:));  
        cost = [cost Inf*ones(Ni,Ni-N)];
		assignment = assignmentoptimal(cost); 	
		assignment(find(assignment==0)) = N+1;
        %perm = sparse(INi(assignment,:));
        perm = sparse(INi(:,assignment)');   
        
        P{k} = perm;
        D2 = [D2 reshape(XF{k}*perm,d2*N,1)];
        %D1((2*i-1):(2*i),:) = XC((2*i-1):(2*i),:)*P_init{i};
        D1 = [D1; XC{k}*perm];
    end    
      	
	leq1 = D1-L-E1;
	leq2 = L-L*C;
	leq3 = D2-M-E2;	
	leq4 = affine*(ones(1,N)*C-ones(1,N));
	
	conv.objVal(iter) = 0.5*norm(C,'fro')^2+lambda1*norm(E1,1)+lambda2*norm_nuc(M)+lambda3*norm(E2,1);
    conv.priRes1(iter) = norm(leq1,'fro');
    conv.priRes2(iter) = norm(leq2,'fro');
	conv.priRes3(iter) = norm(leq3,'fro');
	conv.priRes4(iter) = norm(leq4,'fro');
	conv.priRes(iter) = max( max(conv.priRes1(iter),conv.priRes2(iter)), ...
		max(conv.priRes3(iter),conv.priRes4(iter)) );
	
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










