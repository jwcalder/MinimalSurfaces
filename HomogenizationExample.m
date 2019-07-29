close all
m = 32;
n = 4*m;

[X,Y] = meshgrid(0:1/(n-1):1);
dx = 1/(n-1);
[Xm,Ym] = meshgrid(0:1/(m-1):1);
b = double(rand(m,m) > 0.5);
ae = b + (1-b)*9;
ae = interp2(Xm,Ym,ae,X,Y,'nearest');
f = ones(n,n);
[v,ob1,ob2,ui] = obstacle(1,n);

imagesc(ae);
pause
figure;

disp('PDE Acceleration: Matlab');
u2=LinObs_PDE(ae,f,ob1,ob2,ui,1e6,dx^2,6*pi);
tic;u2=LinObs_PDE_mex(ae,f,ob1,ob2,ui,1e6,dx^2,6*pi);toc; %Mex code is faster

