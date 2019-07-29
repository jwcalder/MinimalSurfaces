close all;
n = 64;
dx = 1/(n-1);
[v,ob1,ob2,ui] = obstacle(1,n);
tol = dx*max(abs(ob1(:)));

disp('PDE Acceleration');
u1 = NonLinObs_PDE(v,ob1,ob2,ui,1e6,tol);
tic;u2 = NonLinObs_PDE_mex(v,ob1,ob2,ui,1e6,tol);toc;

disp('Primal Dual');
u3 = NonLinObs_primaldual(v,ob1,ob2,ui,1e6,tol);
tic;u4 = NonLinObs_primaldual_mex(v,ob1,ob2,ui,1e6,tol);toc;

disp('L1 Penalty method');
tic;u5 = NonLinObs_L1penalty(ob1,ob2,ui,zeros(size(ob1)),1e6,tol,100,500);toc;
tic;u6 = NonLinObs_L1penalty_mex(ob1,ob2,ui,zeros(size(ob1)),1e6,tol,100,500);toc;

