% Author Marwan Torki
% The codes are written according to the paper's description,"Putting Local Features on a Manifold" CVPR 2010
% The paper is co-authored with my advisor Prof. Ahmed Elgammal
% Code is as is and there are no warranties.
% Hope it useful to you.
% Some parts of the m files are witten by their corresponding authors like GB
% feature extraction and the manifold visualization code as well.
% Marwan Torki
% Jul., 18, 2010.
%%
% Here most of the parameters can be edited
clear all;
%Put the folder that has the images you want to compute the manifold for
dataset_path='C:\Users\Marwan\Downloads\demo code\demo code\train9\*.png';
% Number of extracted features per image is an input to the GB code
m=80;
% Reduced dimensionality is the dimension of the y's
Rdim=60;
% flag parameter =0 to use raw gaussian weights in the feature kernel
% flag parameter =1 to use soft correspondencs on gaussian weights using
% Scott and Longuet-Higgins algorithm. For more details look at our other
% paper "one-shot multi-set non-rigid feature-spatial matching" CVPR 2010
f=0;

%%
% compute descriptors with their locations
[descr locs train_files_list]=load_trainSetGB(dataset_path,m);
% generate large weight matrix for all features in the trainset
W=generate_W_blocksone(descr,locs,f);
W=(W+W')/2;
W=W./max(W(:));

instcount=length(locs);

y=compute_embedding(W,Rdim);
% Normalize the y's
for i=1:size(y,1);
    y_n(i,:)=y(i,:)/(norm(y(i,:))+1e-10);
end
%%
% the fllowing line for computing the Hausdorff kernel dissimilarity
% measure on the normailzed y's.
H=generateHausdorff_demo(y_n,m,instcount);
% the fllowing line for computing the squared euclidean distance based on
% the image intensities.
Data_images = collect_images_files(train_files_list,128,128);
IntensityDist=dist2(double(Data_images),double(Data_images));
% You can compare the two matrices which clearly shows the unfold in our
% H better than the IntensityDist matrix.
imagesc(H);
figure;
imagesc(IntensityDist);
%% 
% The following simple lines are the laplacian eigenmaps on the H matrix
% using K nearest Neighbors
ss=1*mean(H(:));
WH=exp(-H/(2*ss));
K=6;
[sorted,index] = sort(WH,'descend');
neighborhood = index(2:(1+K),:);
vals=sorted(2:(1+K),:);
vals=vals';
neighborhood=neighborhood';
D=length(WH);
WNN=zeros(D,D);
for ii=1:D
    jj=neighborhood(ii,:);
    WNN(ii,jj)=vals(ii);
end
WNN=(WNN+WNN')/2;

D=diag(sum(WNN,2));
L=D-WNN;
[V,lamda]=eig(L,D);
Yi=[V(:,2) V(:,3)];
Yres=V(:,2:4);
manifold_image_visualise_m(Yi*30, Data_images, ['image' 'Visualise'], 0.08, [128,128], 0, 0, 1);
