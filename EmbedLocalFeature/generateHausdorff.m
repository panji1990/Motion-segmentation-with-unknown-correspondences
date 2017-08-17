function [H,Yres]=generateHausdorff(V_rd_n,m,icount)
Hausdorff1d=zeros(icount,icount);

        h=1;
        for i=1:icount
            i
                 first=V_rd_n((i-1)*m+1:i*(m),:);
             for j=i+1:icount
                second= V_rd_n((j-1)*m+1:j*m,:);
                Gn1=dist2( first,second);
                Hausdorff1d(i,j)=quantile(min(Gn1'),.5);
                Gn1=dist2( second, first   );
                Hausdorff1d(i,j)=max(Hausdorff1d(i,j),quantile(min(Gn1'),.5));
                h=h+1;
            end
        end
        H=Hausdorff1d;
        H=max(H,H');
        [Yres,eigres]=cmdscale(H);
        