%NonLinObs_primaldual.m - Solves the minimal surface obstacle problem
%
%   -div (nabla u/sqrt(1 + |nabla u|^2)) = f   in  U = (0,1)^2
% 
% subject to boundary condition u=g on partial U
% and obstacle constraints
%
%      ob1 <= u
%
% Other argument:
%
%     ui = initial condition (also encodes g, as ui=g on partial U)
%     T = max number of iterations
%     tol = tolerance.
%
%  Uses the L1-penalty method from 
%
%  Tran, Giang, et al. "An L^1 Penalty Method for General Obstacle Problems." 
%  SIAM Journal on Applied Mathematics 75.4 (2015): 1424-1444.
%
% Author: Jeff Calder, 2019.
%
%
%

function [u,i] = NonLinObs_L1penalty(ob1,ob2,ui,bi,T,tol,lambda,mu)

   pl = true;  %Set to false to supress visualization
   %pl = false;
   s = size(ui);
   ECHECK = 10;

   %Set up parameters
   dx = 1/(s(1)-1);
   dy = 1/(s(2)-1);
   a = 4*pi^2;
   r2 = dx/sqrt(6*a);
   r1 = a*r2;

   %Initialize arrays
   u = ui;
   b = bi;
   v = ob1 - u;

   L = 2/dx^2;
   dt = 1/L/4;
   alpha = (sqrt(L) - sqrt(lambda))/(sqrt(L) + sqrt(lambda));

   err = 1;
   i = 0;
   outer = 1;
   while err > tol & i < T
      outer = outer + 1;

      U = u;
      Uprev = u;
      err2 = 1;
      count2 = 0;
      while err2 > tol
         count2 = count2 + 1;
         i = i+1;
         w = U + alpha*(U - Uprev);
         Uprev = U;

         %Gradient computation
         [ux,uy] = grad(w,dx,dy);
         N = 1./sqrt(1 + ux.^2 + uy.^2);
         G = -div(ux.*N,uy.*N,dx,dy); 

         U = w - dt*(G + lambda*(v - ob1 + w + b));
         U(1,:) = ui(1,:); U(s(1),:) = ui(s(1),:); U(:,1) = ui(:,1); U(:,s(2)) = ui(:,s(2));
         err2 = max(max(abs(U(2:s(1)-1,2:s(2)-1)-Uprev(2:s(1)-1,2:s(2)-1))));
      end
      %disp(count2)
      u = U;
      v = Sp(ob1 - u - b, mu/lambda);
      b = b + u + v - ob1;

      N = 1./sqrt(1 + ux.^2 + uy.^2);
      G = div(ux.*N,uy.*N,dx,dy); 
      F = abs(min(max(G,ob1 - u),ob2 - u));
      err = max(max(F(2:s(1)-1,2:s(2)-1)));

      %Display image
      if pl & mod(outer,100)==0
         m = round(s(1)/64);
         [X,Y]=meshgrid(0:dx:1,0:dy:1);
         surf(X(1:m:end,1:m:end),Y(1:m:end,1:m:end),u(1:m:end,1:m:end));
         %zlim([min(ui(:)),max(ui(:))]);
         drawnow
      end
     
   end
   X = sprintf('Number of Iterations = %d',i);
   disp(X);
end

%Shrink operator
function S = Sp(z,c)
%   S = max(z-c,0);
%   S(z < 0) = z(z<0);
   S = max(z-c,0) + min(z,0);
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


