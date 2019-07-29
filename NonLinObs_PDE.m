%NonLinObs_PDE.m - Solves the minimal surface obstacle problem
%
%   -div (nabla u/sqrt(1 + |nabla u|^2)) = f   in  U = (0,1)^2
% 
% with forcing f, subject to boundary condition u=g on partial U
% and obstacle constraints
%
%      ob1 <= u <= ob2
%
% Other argument:
%
%     ui = initial condition (also encodes g, as ui=g on partial)
%     T = max number of iterations
%     eps = tolerance.
%
% The method used is PDE acceleration
%
%J. Calder and A. Yezzi. PDE Acceleration: A convergence rate analysis and applications to obstacle problems. 2018
%
%M. Benyamin, J. Calder, G. Sundaramoorthi, and A. Yezzi. Accelerated PDE’s for efficient solution of regularized inversion problems. 2018
%
% Author: Jeff Calder, 2018.
%
%

function [u,i] = NonLinObs_PDE(f,ob1,ob2,ui,T,eps)

   pl = true;  %Set to false to supress visualization
   s = size(ui);

   dx = 1/(s(1)-1);
   dy = 1/(s(2)-1);
   [X,Y] = meshgrid(0:dx:1,0:dy:1);
   u = ui;
   uprev = u;
   mu = 1e5;
   dt = 0.8*sqrt(1/2)*dx;
   a = 2*pi;
   
   err = 1;
   i = 0;
   while err > eps & i < T

      %Compute gradient
      [ux,uy] = grad(u,dx,dy);
      N = 1./sqrt(1 + ux.^2 + uy.^2);
      G = div(ux.*N,uy.*N,dx,dy) + f; 

      %Compute residual
      F = abs(min(max(G,ob1 - u),ob2 - u));
      err = max(max(F(2:s(1)-1,2:s(2)-1)));

      %update
      t = uprev;
      uprev = u;
      u = ((2+a*dt)*u - t + dt^2*G)/(1 + a*dt);
      u = min(max(u,ob1),ob2);

      %Boundary condition
      u(1,:) = ui(1,:);u(s(1),:) = ui(s(1),:);u(:,1) = ui(:,1);u(:,s(2)) = ui(:,s(2));

      %Display image
      if pl 
         m = round(s(1)/64);
         surf(X(1:m:end,1:m:end),Y(1:m:end,1:m:end),u(1:m:end,1:m:end));
         %zlim([min(ui(:)),max(ui(:))]);
         drawnow
      end
   
      i = i+1;
   end
   u = uprev;
   X = sprintf('Number of Iterations = %d',i);
   disp(X);
end

%Compute gradient
function [p1,p2] = grad(u,dx,dy)

   s = size(u);
   n = s(1); m = s(2);

   p1 = (u([2:n,n],:) - u)/dx;
   p2 = (u(:,[2:n,n]) - u)/dy;

end

%Compute divergence
function u = div(p1,p2,dx,dy)

   s = size(p1);
   n = s(1); m = s(2);

   p1n = padarray(p1,[1 0],'pre');
   p2n = padarray(p2,[0 1],'pre');
   p1(n,:) = 0;
   p2(:,m) = 0;

   u = (p1 - p1n(1:n,:))/dx + (p2 - p2n(:,1:m))/dy;

end

