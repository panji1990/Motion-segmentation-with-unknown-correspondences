function [I, descr, pos]=GB_desc(imagePath,ndescriptors)
I = im2double(imread(imagePath));
sizeI=size(I)
asp=sizeI(1)/sizeI(2);
I=imresize(I,[200 round(200*asp)]);
if (size(I,3)>1),  % convert the image to grayscale
  I = mean(I,3);
end
%I=imresize(I,[80 80]);

fbr = compute_channels_oe_nms(I);
% rs are the radii for sample points
rs =      [0 4 8 16 32 50];

% nthetas are the number of samples at each radii
nthetas = [1 8 8 10 12 12];

% alpha is the rate of increase for blur
alpha = 0.5;

% beta is the base blur amount
beta = 1;

% Number of descriptors to extract per image
%ndescriptors = 300;

% repulsion radius for rejection sampling
rrep = 5;


% Actually extract Geometric Blur descriptors for each image
[descriptors, pos] = get_descriptors(fbr,ndescriptors,rrep,alpha,beta,rs,nthetas);

descriptors = permute(descriptors,[2 3 1]);            
descriptors = reshape(descriptors,[size(descriptors,1)*size(descriptors,2),size(descriptors,3)])';
descr=descriptors;