function [v,ob1,ob2,ui] = obstacle(ob,n)

   [X,Y] = meshgrid(0:1/(n-1):1);
   z = zeros(n,n);
   o = ones(n,n);
   dx = 1/(n-1);

   switch ob
      case 1  %Obstacle 1 from paper
         g = zeros(n,n);
         g(abs(X-0.6) + abs(Y-0.6) < 0.04) = 5;
         g((X-0.6).^2 + (Y-0.25).^2 < 0.001) = 4.5;
         g(X > 0.075 & X < 0.13 & abs(Y - 0.57) < 1/(n-1)) = 4.5;
         g = g/50;

         ob1 = g;
         ob2 = 1e10*o;
         v = z;      
         ui = g;
         ui(1,:)=0;ui(n,:)=0;ui(:,1)=0;ui(:,n)=0;

      case 2 %Obstacle 2 from paper
         d1 = sqrt((X-0.55).^2 + (Y-0.5).^2)/0.3;
         d2 = sqrt((X-0.1).^2 + (Y-0.5).^2)/0.05;
         g = sqrt(max(0,1-d2.^2)) + sqrt(max(0,1-d1.^2)); 

         ob1 = g;
         ob2 = 1e10*o;
         v = z;
         ui = g;
         ui(1,:)=0;ui(n,:)=0;ui(:,1)=0;ui(:,n)=0;

      case 3 %Double obstacle with forcing
         ob2 = 0.2*ones(size(X));
         ob1 = -min(min(min(X,1-X),Y),1-Y);
         g = zeros(size(X));
         g(0 <= X & X <= 1/6) = 6*X(0 <= X & X <= 1/6);
         g(1/6 < X & X <= 1/3) = 2*(1-3*X(1/6 < X & X <= 1/3));
         g(1/3 < X & X <= 1/2) = 6*(X(1/3 < X & X <= 1/2)-1/3);
         g(1/2 < X & X <= 2/3) = 2*(1-3*(X(1/2 < X & X <= 2/3)-1/3));
         g(2/3 < X & X <= 5/6) = 6*(X(2/3 < X & X <= 5/6)-2/3);
         g(5/6 < X & X <= 1) = 2*(1-3*(X(5/6 < X & X <= 1)-2/3));
         S = abs(X-Y) <= 0.1 & X <= 0.3;
         v = zeros(size(X));
         v(S) = 300;
         I = ~S & X <= 1-Y;
         v(I) = -70*exp(Y(I)).*g(I);
         I = ~S & X > 1-Y;
         v(I) = 15*exp(Y(I)).*g(I);
         
         v = v/10;
         ob1 = ob1/10;
         ob2 = ob2/10;
         ui = ob1;

      otherwise  %Scherk surface
         uscherk = log(cos(X)./cos(Y));
         ui = uscherk;
         ui(2:n-1,2:n-1) = 0;
         v = z;
         ob1 = -1e10*o;
         ob2 = 1e10*o;
   end

end
