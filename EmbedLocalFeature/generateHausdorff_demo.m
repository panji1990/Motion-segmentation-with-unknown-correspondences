function H=generateHausdorff_demo(V_rd_n,m,icount)
Hausdorff1d=zeros(icount,icount);
       for i=1:icount
                 first=V_rd_n((i-1)*m+1:i*(m),:);
             for j=i+1:icount
                second= V_rd_n((j-1)*m+1:j*m,:);
                Gn1=dist2( first,second);
                Hausdorff1d(i,j)=quantile(min(Gn1,[],2),.5);
                Hausdorff1d(i,j)=max(Hausdorff1d(i,j),quantile(min(Gn1),.5));
            end
        end
        H=Hausdorff1d;
        H=max(H,H');
        