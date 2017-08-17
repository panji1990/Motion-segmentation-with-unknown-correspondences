function [accuracy,PN] = Precision1(P_true,P)
K = size(P_true,2);

for j=1:K
cor_num = 0;
tol_num = 0;
PN{j} = P{j}'*P_true{j};
for i=1:K
	temp1 = P{i}*PN{j}; %align to the first permutation
	temp = P_true{i}.*temp1;
	cor_num = cor_num + sum(temp(:));
	temp = P_true{i};
	tol_num = tol_num + sum(temp(:));
end

acc(j) = cor_num/tol_num;
end
[accuracy,ind] = max(acc);
PN = PN{ind};