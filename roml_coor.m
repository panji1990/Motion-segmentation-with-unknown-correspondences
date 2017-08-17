function [PM,D,L,E] = roml_coor(D,lambda)
[d, N] = size(D);
F = d/2;
if(nargin<2)
	lambda = 5/sqrt(N);
end

Aleq = sparse(kron(speye(N),ones(1,N)));
bleq = ones(N,1);
Aeq = sparse(kron(ones(1,N),speye(N)));
beq = ones(N,1);
lb = zeros(N*N,1);
option = optimset('LargeScale','on','Display','off');

% initiate parameters for ADMM
rho = 1e-8; 
%rho = 1e-2;
eta = 1.001;
max_rho = 1e10; maxIter = 1e8;
epsilon = 1e-8;

IN = eye(N);
Id = eye(d);
%I2N = eye(2*N);
X = D;
L = zeros(d,N);
E = zeros(d,N);
Y = zeros(d,N);
for i = 1:F
	PM{i} = sparse(IN);
end
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
	nuc_norm = sum(sigma);
	L = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
	% update E
	temp = D-L+Y/rho;
	E = max(0,temp-lambda/rho)+min(0,temp+lambda/rho);
	% update Pk
	for k=1:F
		xk = [X(2*k-1,:);X(2*k,:)];	
		e = [Id(:,2*k-1) Id(:,2*k)];		
		cost = ((Y'-rho*L'-rho*E')*e*xk)';	
		cost_pos = cost+abs(min(cost(:)));
		%cost_pos = round(1e5*cost_pos);
		%tic
		%assignment = assignmentsuboptimal1(cost_pos);
		tic
 		assignment = assignmentoptimal(cost_pos);	
		toc
		
		P = sparse(IN(assignment,:));
% 		%t1 = toc;

% 		tic        
% 		f = cost(:);
% 		if(max(f)==0)
% 			P = sparse(IN);
% 		else
% 			P = linprog(f',Aleq,bleq,Aeq,beq,lb,[],[],option);
% 			P = round(reshape(P,N,N));
% 			P = sparse(P);
% 		end
% 		t2 =toc;
		
		PM{k} = P;
		D((2*k-1):(2*k),:) = X((2*k-1):(2*k),:)*P;		

% 		test = isequal(perm,P);
	end	
	
	leq = D-L-E;
	stopC = max(max(abs(leq)));	
	
	if(iter==1 || mod(iter,50)==0 || stopC<epsilon)
		disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ...
			',rank=' num2str(rank(L,1e-3*norm(L,2))) ',stopALM=' num2str(stopC,'%2.3e')]);		
		disp(['Obj ' num2str(nuc_norm+lambda*norm(E,1)+trace(Y'*leq)+0.5*rho*norm(leq,'fro')^2)]);
	end
	if(stopC<epsilon)
		break;
	else
		Y = Y + rho*leq;		
		rho = min(max_rho,rho*eta);
	end
end
test = 1;