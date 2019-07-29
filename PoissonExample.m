close all
n = 64;
dx = 1/(n-1);
[X,Y] = meshgrid(0:1/(n-1):1);
o = ones(n,n); z = zeros(n,n);
g = sin(2*pi*X.^2) + cos(2*pi*Y.^2);  %Dirichlet condition
tol = dx^2;  %Tolerance

disp('PDE Acceleration: Matlab');
u1=LinObs_PDE(o,z,-1e10*o,1e10*o,g,1e6,tol);

disp('PDE Acceleration: C Code');
tic;u2=LinObs_PDE_mex(o,z,-1e10*o,1e10*o,g,1e6,tol);toc;  %Mex version is faster

disp('Primal Dual: Matlab');
u3=LinObs_primaldual(z,-1e10*o,1e10*o,g,1e6,tol);

disp('Primal Dual: C Code');
tic;u4=LinObs_primaldual_mex(z,-1e10*o,1e10*o,g,1e6,tol);toc;  %Mex version is faster

