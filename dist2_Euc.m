function n2 = dist2_Euc(x, c)

[dimx, ndata] = size(x);
[dimc, ncentres] = size(c);
if dimx ~= dimc
	error('Data dimension does not match dimension of centres')
end

n2 = sum(c.^2)'*ones(1,ndata)+ones(ncentres,1)*sum(x.^2)-2*c'*x;
    
% this modification is added to avoid generating negative numbers because of numerical errors. 

n2(n2<0)=0;
