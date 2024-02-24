function Test_Poisson

M = 100;
N = 100;
dx = 1 / (M-1);
dy = 1 / (N-1);

 for i = 1:M
     for j = 1 :N   
         
         x = (i-1) * dx;
         y = (j-1) * dy;
         index = (i-1)*M + j;
         
         % On the boundary specify boundary conditions
         if (i==1 || i==M ||j==1 ||j==N)
             rhs(index) = sin(x)*sin(y);
         else
             % If not it is the Laplacian of the exact value
             rhs(index) = -2*sin(x)*sin(y);
         end
         
         Exact(index) = sin(x)*sin(y);
         
     end
 end

 rhs = rhs';
 
 A=  Poisson_LHS(M,N,dx,dy);
 
 p = A\rhs;
 
 norm(p - Exact',2)
 
 plot(Exact,'b');
 hold on
 plot(p,'r o');
 
 
 
 