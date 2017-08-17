clear, clc, close all
%addpath('SSC_1.0');

d1 = 30; d2 = 30;
indn = 0.5; %noise on/off
nx = indn*(randn(1,d1+d2));
ny = indn*(randn(1,d1+d2));
k1 = 2;
b1 = 10;
k2= -1;
b2 = 15;

x1 = (-d1/2):(d1/2-1);
y1 = k1.*x1+b1;

x2 = (-d2/2):(d2/2-1);
y2 = k2.*x2+b2;

X = [x1 x2];
Y = [y1 y2];

X = X+nx;
Y = Y+ny;
figure('Color',[1 1 1]);
plot(X,Y,'b.')
grid on
X = [X;Y];
s = [ones(1,d1) 2*ones(1,d2)];


r = 0; affine = true; outlier = false; rho = 0.7;
alpha = 800;
[missrate1,C1,Grps] = SSC(X,r,affine,alpha,outlier,rho,s);

figure('Color',[1 1 1]);    
X = X';   
[~,ind] = min(missrate1);
for i=1:(d1+d2)
    if(Grps(i,ind) == 2)
        plot(X(i,1),X(i,2),'r.');
        hold on;        
    elseif(Grps(i,ind) == 1)
        plot(X(i,1),X(i,2),'b.');
        hold on; 
        elseif(Grps(i,ind) == 3)
        plot(X(i,1),X(i,2),'r*');
        hold on;
    elseif(Grps(i,ind) == 4)
        plot(X(i,1),X(i,2),'y.');
        hold on;
    elseif(Grps(i,ind) == 5)
        plot(X(i,1),X(i,2),'LineStyle','o');
        hold on;    
    end
end
grid on   

    
