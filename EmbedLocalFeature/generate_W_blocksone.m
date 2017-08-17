function W=generate_W_blocksone(descs,locs,flag)
W=[];
im_count=length(descs);
for i=1:im_count
    X=[locs{i} ones(length(locs{i}),1)];
    X1= comp_affine_invariant_coord(X);
    Glocs=dist2(X1,X1);        
    % uncomment the next line if you want to use gaussian kernels without
    % affine invariant transformation
    %Glocs=dist2(locs{i},locs{i});
    sigma=.01*max(Glocs(:));
    Wlocs=exp(-Glocs/(2*sigma));
    M{i}= .5*Wlocs/sum(Wlocs(:));
    row=[];
for j=1:im_count
    if(i==j)
        row=[row,M{i}];
    else if (i<j)
        Gdescs=dist2(descs{i},descs{j});           
        sigma=.1*max(Gdescs(:));
        Wdescs=exp(-Gdescs/(2*sigma));
        if (flag==1)
            [a b c]=svd(Wdescs);
            Pn=a*eye(size(b))*c';
            Pn=Pn.*(Pn>0.3);
            Cpq{i,j}=Pn/sum(Pn(:));
        else
        Cpq{i,j}=Wdescs/sum(Wdescs(:));
        end
        row=[row, Cpq{i,j}];
        else row=[row, Cpq{j,i}'];
        end
    end

end
 W=[W;row];
end