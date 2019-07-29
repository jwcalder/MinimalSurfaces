%LinObs_PDE.m - Solves the Dirichlet Problem
%
%   -div (ae*grad u) = f   in  U = (0,1)^2
% 
% subject to boundary condition u=g on partial U
% and obstacle constraints
%
%      ob1 <= u <= ob2
%
%  ob1 encodes Dirichlet condition
%
%  Parameters:
%
% ui = Initial condition (also encodes Dirichlet condition)
% T = max number of iterations
% eps = tolerance.
%
% The method used is PDE acceleration
%
%J. Calder and A. Yezzi. PDE Acceleration: A convergence rate analysis and applications to obstacle problems. 2018
%
%M. Benyamin, J. Calder, G. Sundaramoorthi, and A. Yezzi. Accelerated PDE’s for efficient solution of regularized inversion problems. 2018
%
% Author: Jeff Calder, 2018.

function [u,i] = LinObs_PDE(ae,f,ob1,ob2,ui,T,eps,a)

   pl = false;  
   pl = true;  %Turn on to show plot (significnatly slower runtime)
   s = size(ui);

   [X,Y] = meshgrid(1:s(1),1:s(2));
   dx = 1/(s(1)-1);
   dy = 1/(s(2)-1);

   afx = (ae([2:s(1),s(1)],:) + ae)/2;
   afy = (ae(:,[2:s(1),s(1)]) + ae)/2;
   abx = (ae([1,1:s(1)-1],:) + ae)/2;
   aby = (ae(:,[1,1:s(1)-1]) + ae)/2;

   dt = 0.8*sqrt(1/2/max(ae(:)))*dx;

   if nargin == 7
      a = 2*pi;
   end

   u = ui;
   uprev = ui;
   
   err = 1;
   i = 0;
   t = 0;
   while err > eps & i < T  | i < 20

      i = i+1;
     
      uxf = DF(u,1,dx);
      uxb = DB(u,1,dx);
      uyf = DF(u,2,dy);
      uyb = DB(u,2,dy);
      
      Potential(i) = (1/2)*sum(uxf(:).^2 + uyf(:).^2)/s(1)/s(2);
      Kinetic(i) = (1/2)*sum((u(:)-uprev(:)).^2/dt^2)/s(1)/s(2);

      G = uxf.*afx/dx - uxb.*abx/dx + uyf.*afy/dy - uyb.*aby/dy + f; 
      F = abs(min(max(G,ob1 - u),ob2 - u));
      err = max(max(abs(F(2:s(1)-1,2:s(2)-1))));
    
      v = ((2+a*dt)*u - uprev + dt^2*G)/(1+a*dt);
      uprev = u;

      %Obstacle
      u = min(max(v,ob1),ob2);

      %Dirichlet condition
      u(1,:) = ui(1,:);
      u(s(1),:) = ui(s(1),:);
      u(:,1) = ui(:,1);
      u(:,s(2)) = ui(:,s(2));

      %Display graph
      if pl
         m = round(s(1)/64);
         surf(X(1:m:end,1:m:end),Y(1:m:end,1:m:end),u(1:m:end,1:m:end));
         %zlim([min(ui(:)),max(ui(:))]);
         drawnow
      end

   end
   u = uprev;
   X = sprintf('Number of Iterations = %d',i);
   disp(X);
   if pl & 0
      figure
      plot(1:i,Kinetic,'LineStyle','-','LineWidth',2);
      hold
      plot(1:i,Potential,'LineStyle','--','LineWidth',2);
      plot(1:i,Kinetic + Potential,'LineStyle',':','LineWidth',2);
      legend('Kinetic','Potential','Total');
      xlabel('Iterations');
      ylabel('Energy');
   end

end


%Compute centered difference in i with boundary reflection
function A = DC(u,i,dx)

    A = (DF(u,i,dx) + DB(u,i,dx))/2;
    
end

%Compute forward difference in i with boundary reflection
function A = DF(u,i,dx)

    s = size(u);
    switch i
        case 1
            A = (u([2:s(1),s(1)-1],:) - u)/dx;
        case 2
            A = (u(:,[2:s(2),s(2)-1]) - u)/dx;
    end
    
end

%Compute backward difference in i with boundary reflection
function A = DB(u,i,dx)

    s = size(u);
    switch i
        case 1
            A = (u - u([2,1:s(1)-1],:))/dx;
        case 2
            A = (u - u(:,[2,1:s(2)-1]))/dx;
    end
    
end

