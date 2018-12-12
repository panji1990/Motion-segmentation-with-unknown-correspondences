% Author: Pan Ji, ANU
% Copyright reserved!
function [sift_acc, P_sift]= sift_matching_hungarian(Fea,P)
F = length(Fea);
N = size(Fea{1},2);
if(nargin<2)
	for i=1:F
		P{i} = speye(N);
	end
end
PN = P{1}';
for i=1:F
	P{i} = P{i}*PN;
end
euc_dist = [];
P_sift = [];
baseline = 1:N;
num = 0;
tol = N*(F-1);
IN = eye(N);

P_sift{1} = speye(N);
for i=1:F-1    
	cost = dist2_Euc(Fea{i+1},Fea{1});	
    euc_dist{i} = cost;	 
	assignment = assignmentoptimal(double(cost));
	baseline_p = baseline*P{i+1};
	num = num+sum(assignment(:)==baseline_p(:));
	P_sift{i+1} = sparse(IN(:,assignment));
end
sift_acc = num/tol;
