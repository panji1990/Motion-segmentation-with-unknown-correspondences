
function Data_images = collect_images_files(train_files,Nx,Ny)

%dataset_path='/Users/ahmed/Documents/Research/Marwan Code/work/Localization Expr/train2/*.jpg';

%train_files=fuf(dataset_path,'detail');

Data_images=uint8(zeros(size(train_files,1),Nx*Ny));

for i=1:length(train_files)
   im=imresize(imread(train_files{i}),[Nx,Ny]);
   if size(im,3)==3,
    im=rgb2ycbcr(im);
   end
   
   Data_images(i,:)=reshape(im(:,:,1),1,Nx*Ny);
    
end

%save shape_train_images Data_images
