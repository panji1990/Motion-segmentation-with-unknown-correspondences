function [descr locs,train_files]=load_trainSetGB(dataset_path,m)
train_files=fuf(dataset_path,'detail');
for i=1:length(train_files)
     i
     [I descr{i} locs{i}]=GB_desc(train_files{i},m);
 end