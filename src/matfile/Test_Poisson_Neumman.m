function Test_Poisson

clear all
M = 120;
N = 120;
dx = pi / (M-2);
dy = pi / (N-2);

 for i = 1:M
     for j = 1 :N   
         
         x = (i-1 - 0.5) * dx;
         y = (j-1 - 0.5) * dy;
         
         index = glidx(i,j,M,N);
         
         % On the boundary specify boundary conditions
         if (i==1 || i==M ||j==1 ||j==N)
             rhs(index) = 0;
         else
             % If not it is the Laplacian of the exact value
             rhs(index) = -2*cos(x) *cos(y);
         end
         
         % We don't solve it
         if (i==2 && j==2)
             rhs(index) = 0;%cos(x) * cos(y);
         end

         Exact(index) = cos(x) * cos(y);
         
         Exact_Phi(i,j) = cos(x) * cos(y);
     end
 end
 
 Ucont(1:M,1:N) =0;
 P_RHS = Poisson_RHS_Neumann(Ucont,Ucont,dx,dy,0);
 A = Poisson_LHS_Neumann(M,N,dx,dy); 
 phi = A\rhs';
 
 % difference only by a constant
 index = glidx(2,2,M,N);
 phi = phi + Exact(index);
 
 figure(1)
 plot(Exact,'b');
 hold on
 plot(phi,'r o')
 
  %% At the face center
 
 
 for i = 1:M-1
     for j= 1:N         
         
         x = (i-1) * dx;
         y = (j-1 - 0.5) * dy;
         
         G_x(i,j) = -sin(x)*cos(y);
         
     end
 end

 
 
 for i = 1:M
     for j= 1:N-1         
 
         x = (i-1 - 0.5) * dx;
         y = (j-1) * dy;
 
         G_y(i,j) = -cos(x)*sin(y);
         
     end
 end
 
 [U_x_new U_y_new P_new] = Update_Solution(Ucont, Ucont, Ucont, phi,dx,dy,0.1,0);
 
 error(1:M,1:N)  = 0;
 for i = 2:M-1
     for j = 2:N-1
         
         error(i,j) = P_new(i,j) - Exact_Phi(i,j); 
     end
 end
 
 figure(2)
 mesh(error);