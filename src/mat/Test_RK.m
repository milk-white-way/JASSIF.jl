function Runge_Kutta()

pseudot = 0;
alpha = [1/4; 1/3; 1/2; 1];

% Form BCS
Ucont_x = 1; 
dt = 0.01;
% Assign first guess
U_p_x = Ucont_x;
U_im_x = Ucont_x;

tol = 1e-8;

e = norm(Ucont_x);

 while (e > tol && pseudot < 15 )     
          
     U_im_x = U_p_x;
     
      for stage=1:4
          
          RHS_x = (U_im_x - 2 ) * (U_im_x - 3);
          U_im_x = U_p_x + alpha(stage) * dt * RHS_x;
          
      end     
     
     e = norm(U_p_x - U_im_x,inf); 
    
     U_p_x = U_im_x;
     
     fprintf('subitr = %d, Convergence of momentum = %8.6f \n',pseudot,e);
     
     pseudot = pseudot+1;    
          
 end 
 U_p_x
