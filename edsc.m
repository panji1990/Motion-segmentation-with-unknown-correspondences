function [Missrate, C, grp,time] = edsc(X,s,lambda,affine,outlier,Dim,alpha)
if(nargin<5)
	outlier = false;
end
if(nargin<4)
	affine = false;
end
if(nargin<3)
	lambda = 100;
	lambda1 = 100;
	lambda2 = 50;
elseif(length(lambda)==2)
	lambda1 = lambda(1);
	lambda2 = lambda(2);
elseif(length(lambda)==1)
	lambda1 = lambda;
	lambda2 = lambda/2;
end
time = 0;
[d,N] = size(X);
K = max(s);

if(~outlier)
	lxtx = lambda*(X'*X);
	if(~affine) %linear
		C = (eye(N)+lxtx)\lxtx;
	else %affine
% 		tic
% 		B = (eye(N)+lxtx)\lxtx;
% 		A = (eye(N)+lxtx)\ones(N,1);
% 		lhand = (A'*A) + lambda*((X*A)'*(X*A)) - 2*(ones(1,N)*A);
% 		rhand = B'*A + lambda*(X'-B'*X')*X*A - (B'*ones(N,1)-ones(N,1));
% 		Y = rhand/lhand;
% 		C = B - A*Y';
% 		t1 = toc;
		
		%test	
		tic
		One = ones(1,N);
		V0 = eye(N);
		V = null(One);
		lhand = V'*V+V'*lxtx*V;
		rhand = -V'*V0;
		phi = lhand\rhand;
		C = V0+V*phi;	
		time = toc;
	end
else
	lxtx = lambda1*(X'*X);
	if(~affine) %linear
		% ALM
		% parameters
		epsilon = 1e-8;
		rho = 1e-6;
		%rho = 0.2;
		maxIter = 1e4;
		eta = 3;
		max_rho = 1e10;
		% initialize
		C = zeros(N,N);
		E = zeros(d,N);
		Y = zeros(d,N);
		
		iter = 0;
		while(iter<maxIter)
			iter = iter+1;
			% update C
			C = (eye(N)+lxtx+rho*(X'*X))\(lxtx+X'*Y+rho*X'*(X-E));
			% update E
			xmxc = X-X*C;
			temp = xmxc+Y/rho;
			E = max(0,temp - lambda2/rho)+min(0,temp + lambda2/rho);
			
			leq = xmxc - E;
			stpC = max(max(abs(leq)));
			if(iter == 1 || mod(iter,50)==0 || stpC<epsilon)
				disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ',stopALM=' num2str(stpC,'%2.3e')]);
			end
			if(stpC<epsilon)
				break;
			else
				Y = Y + rho*leq;
				rho = min(max_rho,rho*eta);
			end
		end
	else %affine
		% ALM
		% parameters
		epsilon = 1e-8;
		rho = 1e-6;
		maxIter = 1e4;
		eta = 1.1;
		max_rho = 1e10;		
		% initialize
		C = zeros(N,N);
		E = zeros(d,N);
		Y1 = zeros(d,N);
		Y2 = zeros(N,1);
		
		iter = 0;
		while(iter<maxIter)
			iter = iter+1;
			% update C
			C = (eye(N)+lxtx+rho*(X'*X+ones(N,N)))\(lxtx+X'*Y1-ones(N,1)*Y2'+rho*(X'*(X-E)+ones(N,N)));
			% update E
			xmxc = X-X*C;
			temp = xmxc+Y1/rho;
			E = max(0,temp - lambda2/rho)+min(0,temp + lambda2/rho);
			
			leq1 = xmxc - E;
			leq2 = C'*ones(N,1)-ones(N,1);
			stpC = max( max(max(abs(leq1))),max(max(abs(leq2))) );
			if(iter == 1 || mod(iter,50)==0 || stpC<epsilon)
				disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ',stopALM=' num2str(stpC,'%2.3e')]);
			end
			if(stpC<epsilon)
				break;
			else
				Y1 = Y1 + rho*leq1;
				Y2 = Y2 + rho*leq2;
				rho = min(max_rho,rho*eta);
			end
		end
	end
end

% W = abs(C)+abs(C');
% grp = ncutW(W,K);
[grp,t2] = post_proC(C,K,Dim,alpha);
%time = t1;
Missrate = ErrorRate(grp,s);










