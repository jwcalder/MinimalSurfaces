%NonLinObs_primaldual.m - Solves the minimal surface obstacle problem
%
%   -div (nabla u/sqrt(1 + |nabla u|^2)) = f   in  U = (0,1)^2
% 
% subject to boundary condition u=g on partial U
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
% The method used is an improved version of the primal dual method from
%
%     Zosso, Dominique, et al. "An Efficient Primal-Dual Method for the Obstacle Problem." 
%     Journal of Scientific Computing 73.1 (2017): 416-437.
%
% Author: Jeff Calder, 2018.
%
%
%

function u = NonLinObs_primaldual(f,ob1,ob2,ui,T,eps)

   pl = true;  %Set to false to supress visualization
   s = size(ui);
   ECHECK = 10;

   %Set up parameters
   dx = 1/(s(1)-1);
   dy = 1/(s(2)-1);
   NB = round(log(eps*dx^2)/log(1/2)); %Number of bisection iterations
   a = 4*pi^2;
   r2 = dx/sqrt(6*a);
   r1 = a*r2;

   %Initialize arrays
   u = ui;
   uprev = u;
   ubar = u;
   ubarprev = u;
   p1 = zeros(s);
   p2 = zeros(s);

   err = 1;
   i = 0;
   while err > eps & i < T

      i = i+1;
      %Dual update
      [ux,uy] = grad(ubar,dx,dy);
      q1 = p1 + r1*ux;
      q2 = p2 + r1*uy;
      N = sqrt(q1.^2 + q2.^2 + 1e-10); 
      q1 = q1./N;
      q2 = q2./N;
      
      %Bisection search
      alpha1 = zeros(s);
      alpha2 = min(ones(s),N);
      alpha = (alpha1 + alpha2)/2;
      for j=1:NB
         g = r1^2*alpha.^2 - (1-alpha.^2).*(alpha - N).^2;
         I = double(g>0);
         alpha1 = I.*alpha1 + (1-I).*alpha;
         alpha2 = I.*alpha  + (1-I).*alpha2;
         alpha = (alpha1 + alpha2)/2; 
      end
      p1 = alpha.*q1;
      p2 = alpha.*q2;
      
      %Primal update
      u = u + r2*div(p1,p2,dx,dy) + r2*f;

      %Obstacle and boundary conditions
      u = min(max(u,ob1),ob2);
      u(1,:) = ui(1,:);u(s(1),:) = ui(s(1),:);u(:,1) = ui(:,1);u(:,s(2)) = ui(:,s(2));

      %Error computation
      if mod(i,ECHECK) == 0
         N = 1./sqrt(1 + ux.^2 + uy.^2);
         G = div(ux.*N,uy.*N,dx,dy) + f; 
         F = abs(min(max(G,ob1 - u),ob2 - u));
         err = max(max(F(2:s(1)-1,2:s(2)-1)));
      end

      ubarprev = ubar;
      ubar = 2*u - uprev;
      uprev = u;

      %Display image
      if pl
         m = round(s(1)/64);
         [X,Y]=meshgrid(0:dx:1,0:dy:1);
         surf(X(1:m:end,1:m:end),Y(1:m:end,1:m:end),ubar(1:m:end,1:m:end));
         %zlim([min(ui(:)),max(ui(:))]);
         drawnow
      end
     
   end
   u = ubarprev;
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


