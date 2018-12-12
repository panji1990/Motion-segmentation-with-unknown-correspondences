% Author: Pan Ji, ANU
% Copyright reserved!
function [XP,P] = permutated(X,option)
if(nargin<2)
	option = 'Point'; % Feature
end

if(strcmpi(option,'Point'))
	[d,N] = size(X);
	F = d/2;
	
	IN = eye(N);
	for k = 1:F
		perm = randperm(N);
		Pf = sparse(IN(perm,:));
		XP(2*k-1,:) = X(2*k-1,:)*Pf;
		XP(2*k,:) = X(2*k,:)*Pf;
		P{k} = Pf';
	end
% 	for k=1:3
%  	P{1} = speye(N);
%   	XP(1:2,:) = X(1:2,:);	
% 	end
	
elseif(strcmpi(option,'Feature'))
	F = size(X,2);
	N = size(X{1},2);
	IN = eye(N);
	%XP{1} = double(X{1});
	%P{1} = sparse(IN);
% 	for k=1:floor(F/2)
% 		P{k} = sparse(IN);
% 		XP{k} = double(X{k});
% 	end
	for k = 1:F
		perm = randperm(N);
		Pf = sparse(IN(perm,:));
		XP{k} = double(X{k})*Pf;
		P{k} = Pf';
	end
else
	disp('Option: Point or Feature!');
	return;
end
