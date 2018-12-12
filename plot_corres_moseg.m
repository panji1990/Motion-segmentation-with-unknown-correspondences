function [XC,XF] = plot_corres_moseg(video,locs,descrs,PM,grp)
F = length(locs); %the number of frames
N = size(PM{1},2);  %the number of inlier points
n = size(grp,2); %the number of motions

XC = [];
XF = {};
for f=1:F
	XC(2*f-1:2*f,:) = double(locs{f})*PM{f};
	%XF{f} = double(descrs{f})*PM{f};
end

%figure(100)
for f=1:F
%f_show = floor(F/2);
figure
f_show = f;
imshow(uint8(video{f_show}))
for i=1:N
	hold on
    if(grp(i,1))
		plot(XC(2*f_show-1,i),XC(2*f_show,i),'rx','MarkerSize',15,'LineWidth',3)
	else
		plot(XC(2*f_show-1,i),XC(2*f_show,i),'bo','MarkerSize',15,'LineWidth',3)
	end
end
end
figure(101)
%for i=10:N
for f=1:F
	imshow(uint8(video{f}))
	hold on
	plot(XC(2*f-1,1),XC(2*f,1),'rx','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,2),XC(2*f,2),'ro','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,3),XC(2*f,3),'r*','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,4),XC(2*f,4),'r+','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,5),XC(2*f,5),'rd','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,6),XC(2*f,6),'rs','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,7),XC(2*f,7),'r^','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,8),XC(2*f,8),'rv','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,9),XC(2*f,9),'rp','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,10),XC(2*f,10),'r<','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,11),XC(2*f,11),'bx','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,12),XC(2*f,12),'bo','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,13),XC(2*f,13),'b*','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,14),XC(2*f,14),'b+','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,15),XC(2*f,15),'bd','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,16),XC(2*f,16),'bs','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,17),XC(2*f,17),'b^','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,18),XC(2*f,18),'bv','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,19),XC(2*f,19),'bp','MarkerSize',10,'LineWidth',2)
	plot(XC(2*f-1,20),XC(2*f,20),'b<','MarkerSize',10,'LineWidth',2)
% 	plot(XC(2*f-1,21),XC(2*f,21),'ms','MarkerSize',10,'LineWidth',2)
% 	plot(XC(2*f-1,22),XC(2*f,22),'m^','MarkerSize',10,'LineWidth',2)
% 	plot(XC(2*f-1,23),XC(2*f,23),'mv','MarkerSize',10,'LineWidth',2)
% 	plot(XC(2*f-1,24),XC(2*f,24),'mp','MarkerSize',10,'LineWidth',2)
% 	plot(XC(2*f-1,25),XC(2*f,25),'m<','MarkerSize',10,'LineWidth',2)
	%pause(0.01)
end
%end