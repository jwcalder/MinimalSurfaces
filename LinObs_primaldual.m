%LinObs_primaldual_mex.c - Solves the Dirichlet Problem
%
%   -Delta u = f   in  U = (0,1)^2
% 
% subject to boundary condition u=g on partial U
% and obstacle constraints
%
%      ob1 <= u <= ob2
%
% ui = initial condition (also encodes g, as ui=g on partial)
% T = max number of iterations
% eps = tolerance.
%
% The method used is the Primal Dual method from
%
% Zosso, Dominique, et al. "An Efficient Primal-Dual Method for the Obstacle Problem." 
% Journal of Scientific Computing 73.1 (2017): 416-437.
%
% Author: Jeff Calder, 2018.

function u = LinObs_primaldual(f,ob1,ob2,ui,T,eps)

   pl = false;  
   pl = true;  %Turn on to show plot (significnatly slower runtime)
   s = size(ui);

   [X,Y] = meshgrid(1:s(1),1:s(2));
   dx = 1/(s(1)-1);
   dy = 1/(s(2)-1);
   u = ui;
   ubar = u;
   p1 = zeros(s);
   p2 = zeros(s);
   a = 4*pi^2;
   r2 = dx/sqrt(6*a);
   r1 = a*r2;

   err = 1;
   i = 0;
   while err > eps & i < T | i < 20

      %Dual update
      [q1,q2] = grad(ubar);
      p1 = (p1 + r1*q1/dx)/(1+r1);
      p2 = (p2 + r1*q2/dy)/(1+r1);
     
      uprev = u;
      u = u + r2*(div(p1/dx,p2/dy) + f);
      u = min(max(u,ob1),ob2);
      u(1,:) = ui(1,:);
      u(s(1),:) = ui(s(1),:);
      u(:,1) = ui(:,1);
      u(:,s(2)) = ui(:,s(2));

      ubarprev = ubar;
      ubar = 2*u - uprev;

      E = div(q1/dx^2,q2/dy^2) + f;
      F = abs(min(max(E,ob1 - u),ob2 - u));
      err = max(max(abs(F(2:s(1)-1,2:s(2)-1))));

      %Display image
      if pl
         m = round(s(1)/64);
         surf(X(1:m:end,1:m:end),Y(1:m:end,1:m:end),u(1:m:end,1:m:end));
         zlim([min(ui(:)),max(ui(:))]);
         drawnow
      end
    
      i = i+1;
   end
   u = ubarprev;
   X = sprintf('Number of Iterations = %d',i);
   disp(X);
end


%Compute gradient
function [p1,p2] = grad(u)

   s = size(u);
   n = s(1); m = s(2);

   p1 = [diff(u);zeros(1,m)];
   p2 = [diff(u,1,2),zeros(n,1)];

end

%Compute divergence
function u = div(p1,p2)

   s = size(p1);
   n = s(1); m = s(2);

   a = diff(p1);
   a(n-1,:)=[];
   a = [p1(1,:);a;-p1(n-1,:)];

   b = diff(p2,1,2);
   b(:,m-1)=[];
   b = [p2(:,1),b,-p2(:,m-1)];

   u = a + b;

end
