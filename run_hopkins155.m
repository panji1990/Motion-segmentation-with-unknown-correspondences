clear, close all

addpath(genpath(pwd));
cd('HopkinsOutdoor');
warning off
file = dir;
flag = 0;
Rdim = 60;
for i = 1:(length(file)-2)
	close all
	filepath = file(i+2).name
	eval(['cd ' filepath]);
	vname = [filepath '.avi'];
	mname = [filepath '_truth.mat'];
	load(mname);
	n = max(s); %number of clusters
	N = size(x,2); %number of points
	F = size(x,3); % number of frames
	D = 2*F;
	Y = reshape(permute(y(1:2,:,:),[1 3 2]),D,N);
 	X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N); % normalized x
	if(isempty(strfind(vname,'kanatani')))
		vObj = VideoReader(vname);	
		nFrames = vObj.NumberOfFrames;
	end	
	Fea = [];descr={};locs={}; frames = {};
		
	for k = 1:F
		if(isempty(strfind(vname,'kanatani')))
			img = read(vObj,k);
		else
			video = mmread(vname,k);
			img = video.frames.cdata;
		end
		if(numel(size(img))>2) % is rgb or not
			img = rgb2gray(img);
		end
		frames{k} = img;
		%img = img(1:2:end,1:2:end);
		
		figure(1)
		%figure('Color',[1 1 1])
		imshow(img)
		imwrite(img,[filepath '_' sprintf('%02d',k) '.jpg'],'jpg')
		
		
		if(~isempty(strfind(filepath,'cars2_06'))||~isempty(strfind(filepath,'cars2_07'))...
				||~isempty(strfind(filepath,'1RT2TCRT_B'))||~isempty(strfind(filepath,'kanatani')))
 			f = Y(2*k-1:2*k,:);
		else
			f = [Y(2*k-1,:);height-Y(2*k,:)];
		end

		hold on
		plot(f(1,:),f(2,:),'rx')
		
		[test,d,inf] = vl_covdet(single(img),'frames',f,'descriptor','SIFT','verbose');%'EstimateAffineShape', true,'EstimateOrientation', true,,'PeakThreshold',30 ,'verbose'	
		h3 = vl_plotsiftdescriptor(d,test) ;
        set(h3,'color','g') ;
		Fea{k}   = double(d);
		descr{k}   = double(d)';
		locs{k} = X(2*k-1:2*k,:)';
	end
		
		
  	cd ..;
	[XFP,P] = permutated(Fea,'Feature');	
	
	%% ------Embedded Features-----------
	W=generate_W_blocksone(descr,locs,flag);
	W=(W+W')/2;
	W=W./max(W(:));	
	instcount=length(locs);	
	fea_embed=compute_embedding(W,Rdim);
	% Normalize the y's
	fea_embed_n = normr(fea_embed);
	for k=1:F
		XF_embed{k} = fea_embed_n(N*(k-1)+1:N*k,:)';
	end
	
	XCP = [];
	XFP_embed=[];
	for j=1:F
		XCP(2*j-1:2*j,:) = X(2*j-1:2*j,:)*P{j}'; % permutate the coordinate with the same P
		XFP_embed{j} = XF_embed{j}*P{j}'; % permutate the embeded features with the same P
	end
	
	%% -----sift matching using Hungarian algorithm----	
	[Acc_sift(i),P_sift] = sift_matching_hungarian(XFP,P);
	XC_sift = [];
	for k=1:F
		XC_sift(2*k-1:2*k,:) = XCP(2*k-1:2*k,:)*P_sift{k};
	end
	[~,~,grp_sift] = edsc(XC_sift,s,5,true,false,4,4); % SIFT Matching+EDSC
	[sift_acc,PN_sift] = Precision1(P,P_sift);
	grp_sift=(grp_sift'*PN_sift)';
	Missrate_sift_edsc(i) = ErrorRate(grp_sift,s);
	
	[~,~,~,grp_sift_lrr] = solve_lrr(XC_sift,s,4,false,false,0,true);
	grp_sift_lrr=(grp_sift_lrr'*PN_sift)';
	Missrate_sift_lrr(i) = ErrorRate(grp_sift_lrr,s);
	
	[~,~,~,grp_sift_ssc] = SSC(XC_sift,0,true,800,false,0.7,s);
	grp_sift_ssc=(grp_sift_ssc'*PN_sift)';
	Missrate_sift_ssc(i) = ErrorRate(grp_sift_ssc,s);
	
	%% -----ROML with SIFT features
	lambda = 5/sqrt(N);
	[PM_roml_sift,XM,E] = roml_feature(XFP,N,lambda,P_sift);%	
	[Acc_roml_sift(i),PN_roml_sift] = Precision1(P,PM_roml_sift);
	XC_roml_sift = [];
	for k=1:F
		XC_roml_sift(2*k-1:2*k,:) = XCP(2*k-1:2*k,:)*PM_roml_sift{k};
	end
	[~,~,grp_roml_sift] = edsc(XC_roml_sift,s,5,true,false,4,4); % ROML+EDSC
	grp_roml_sift=(grp_roml_sift'*PN_roml_sift)';
	Missrate_roml_sift_edsc(i) = ErrorRate(grp_roml_sift,s);
	
	[~,~,~,grp_roml_sift_lrr] = solve_lrr(XC_roml_sift,s,4,false,true,0,true);
	grp_roml_sift_lrr=(grp_roml_sift_lrr'*PN_roml_sift)';
	Missrate_roml_sift_lrr(i) = ErrorRate(grp_roml_sift_lrr,s);
	
	[~,~,~,grp_roml_sift_ssc] = SSC(XC_roml_sift,0,true,800,false,0.7,s);
	grp_roml_sift_ssc=(grp_roml_sift_ssc'*PN_roml_sift)';
	Missrate_roml_sift_ssc(i) = ErrorRate(grp_roml_sift_ssc,s);
	
	%% ----- ROML with coordinates
% % 	lambda = 5/sqrt(N);
% % 	[PM_roml_coor,~] = roml_coor(XCP,lambda);	
% % 	[Acc_roml_coor(i),PN_roml_coor] = Precision1(P,PM_roml_coor);
% % 	XC_roml_coor = [];
% % 	for k=1:F
% % 		XC_roml_coor(2*k-1:2*k,:) = XCP(2*k-1:2*k,:)*PM_roml_coor{k};
% % 	end
% % 	[~,~,grp_roml_coor] = edsc(XC_roml_coor,s,120,true,false,4,4); % ROML+EDSC
% % 	grp_roml_coor=(grp_roml_coor'*PN_roml_coor)';
% % 	Missrate_roml_coor(i) = ErrorRate(grp_roml_coor,s);
% 	
	%% -----ROML with Embedded features----
	lambda = 1/sqrt(N);
	[PM_roml_embed,XM,E] = roml_feature(XFP_embed,N,lambda,P_sift);%	
	[Acc_roml_embed(i),PN_roml_embed] = Precision1(P,PM_roml_embed);
	XC_roml_embed = [];
	for k=1:F
		XC_roml_embed(2*k-1:2*k,:) = XCP(2*k-1:2*k,:)*PM_roml_embed{k};
	end
	[~,~,grp_roml_embed] = edsc(XC_roml_embed,s,5,true,false,4,4); % ROML+EDSC
	grp_roml_embed=(grp_roml_embed'*PN_roml_embed)';
	Missrate_roml_embed_edsc(i) = ErrorRate(grp_roml_embed,s);
	
	[~,~,~,grp_roml_embed_lrr] = solve_lrr(XC_roml_embed,s,4,false,true,0,true);
	grp_roml_embed_lrr=(grp_roml_embed_lrr'*PN_roml_embed)';
	Missrate_roml_embed_lrr(i) = ErrorRate(grp_roml_embed_lrr,s);
	
	[~,~,~,grp_roml_embed_ssc] = SSC(XC_roml_embed,0,true,800,false,0.7,s);
	grp_roml_embed_ssc=(grp_roml_embed_ssc'*PN_roml_embed)';
	Missrate_roml_embed_ssc(i) = ErrorRate(grp_roml_embed_ssc,s);
	
	%% -----JEPS--------
	lambda1 = 0.05;
	lambda2 = 1;
	lambda3 = lambda2*5/sqrt(N);
	[C,PM,L,M,D1,D2,E1,E2] = JEPS_feature(XCP,XFP,lambda1,lambda2,lambda3,true,P_sift);%	
	[Acc_jeps(i),PN] = Precision1(P,PM);	
	grp_jeps = post_proC(C,max(s));
	grp_jeps = (grp_jeps'*PN)';
	Missrate_jeps(i) = ErrorRate(grp_jeps,s);
	
	for k=1:F
		XC_jeps(2*k-1:2*k,:) = XCP(2*k-1:2*k,:)*PM{k};
	end
	
	[~,~,grp_JEPS] = edsc(XC_jeps,s,5,true,false,4,4); % JEPS+EDSC
	grp_JEPS=(grp_JEPS'*PN)';
	Missrate_jeps_edsc(i) = ErrorRate(grp_JEPS,s);
	
	[~,~,~,grp] = solve_lrr(XC_jeps,s,4,false,true,0,true);
	grp=(grp'*PN)';
	Missrate_jeps_lrr(i) = ErrorRate(grp,s);
	
	[~,~,~,grp] = SSC(XC_jeps,0,true,800,false,0.7,s);
	grp=(grp'*PN)';
	Missrate_jeps_ssc(i) = ErrorRate(grp,s);
	
	% ----EDSC-----
	Missrate_edsc(i) = edsc(X,s,5,true,false,4,4);
	Missrate_lrr(i) = solve_lrr(X,s,4,false,false,0,true);
 	Missrate_ssc(i) = SSC(X,0,true,800,false,0.7,s);

    save('D:\Reference_Code\Subspace Clustering Algorithms\Clustering_without_matching\output\result.mat','Acc_sift','Missrate_sift_edsc','Missrate_sift_lrr','Missrate_sift_ssc',...
		'Acc_roml_sift','Missrate_roml_sift_edsc','Missrate_roml_sift_lrr','Missrate_roml_sift_ssc',...
 		'Acc_roml_embed','Missrate_roml_embed_edsc','Missrate_roml_embed_lrr','Missrate_roml_embed_ssc',...
		'Acc_jeps','Missrate_jeps','Missrate_edsc','Missrate_lrr','Missrate_ssc');
	save('output\cars2_07_g12_result.mat','grp_sift_ssc','grp_roml_sift_ssc','grp_roml_embed_ssc','grp_jeps','Y','frames',...
		'P','PM_sift','PN_sift','PM_roml_sift','PN_roml_sift','PM_roml_embed','PN_roml_embed','PM_jeps','PN_jeps');

end

avg_sift = mean(Acc_sift);
avg_miss_sift_edsc = mean(Missrate_sift_edsc);
avg_miss_sift_lrr = mean(Missrate_sift_lrr);
avg_miss_sift_ssc = mean(Missrate_sift_ssc);
 
avg_roml_embed = mean(Acc_roml_embed);
avg_miss_roml_embed_edsc = mean(Missrate_roml_embed_edsc);
avg_miss_roml_embed_lrr = mean(Missrate_roml_embed_lrr);
avg_miss_roml_embed_ssc = mean(Missrate_roml_embed_ssc);
  
avg_roml_sift = mean(Acc_roml_sift);
avg_miss_roml_sift_edsc = mean(Missrate_roml_sift_edsc);
avg_miss_roml_sift_lrr = mean(Missrate_roml_sift_lrr);
avg_miss_roml_sift_ssc = mean(Missrate_roml_sift_ssc);

avg_jeps = mean(Acc_jeps);
avg_miss_jeps = mean(Missrate_jeps);

avg_miss_edsc = mean(Missrate_edsc);
avg_miss_lrr = mean(Missrate_lrr);
avg_miss_ssc = mean(Missrate_ssc);

save('D:\Reference_Code\Subspace Clustering Algorithms\Clustering_without_matching\output\mean_result.mat','avg_sift','avg_miss_sift_edsc','avg_miss_sift_lrr','avg_miss_sift_ssc',...
	'avg_roml_embed','avg_miss_roml_embed_edsc','avg_roml_sift','avg_miss_roml_sift_edsc','avg_miss_roml_sift_lrr','avg_miss_roml_sift_ssc',...
	'avg_jeps','avg_miss_jeps','avg_miss_edsc','avg_miss_lrr','avg_miss_ssc');
	

