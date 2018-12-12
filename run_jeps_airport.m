% Author: Pan Ji, ANU
% Copyright reserved!
clear, close all
%cd('D:\Reference_Code\Subspace Clustering Algorithms\Clustering_without_matching')
addpath(pwd);
addpath('common');
cd('airportdataset\052');

F = 40; % the number of frames
locs = {}; XC = {};
descrs = {}; XF = {};
video = {};
for f=1:F
	name = ['img-000000' sprintf('%02d',f) '.jpg'];
    img = rgb2gray(imread(name));
	video{f} = img;
	figure(1)
	imshow(img)
    [loc,descr] = vl_sift(single(img),'PeakThresh',5,'EdgeThresh',5);
	hold on
	plot(loc(1,:),loc(2,:),'rx')
	%h3 = vl_plotsiftdescriptor(descr,loc) ;
	%set(h3,'color','g') ;
	locs{f} = double(loc(1:2,:));
	XC{f} = double(loc(1:2,:))/sqrt(1920^2+1080^2);
	descrs{f} = double(descr);
	XF{f} = normc(double(descr));
end

for f=1:F-1
	[matches,scores] = vl_ubcmatch(descrs{1},descrs{f+1});
    num_match(f) = size(matches,2);
end
N = min(num_match)
param.lambda1 = 0.4; param.lambda2 = 1; param.lambda3 = 2/sqrt(N); affine = true;
param.rho = 1e-6;
param.eta = 1.01;
param.epsilon = 1e-6;
tic
[C,PM,L,M,D1,D2,E1,E2,conv] = JEPS_Hungarian(XC,XF,N,affine,param);
time = toc;
grp = post_proC(C,2);

[XCP,XFP] = plot_corres_moseg(video,locs,descrs,PM,grp);
