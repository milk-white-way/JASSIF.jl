function [U_im_x U_im_y] =  Momentum_Solver()

M = 30;
N =30;

myzero = rand(M,N);
 % Set initial guess
 x0.dU_x    = myzero;
 x0.dU_y    = myzero;
 x0.Ucont_x = myzero;
 x0.Ucont_y = myzero;
 x0.Ucat_x  = myzero;
 x0.Ucat_y  = myzero;
 x0.Ubcs_x  = myzero;
 x0.Ubcs_y  = myzero;
 x0.Pressure= myzero;
 x0.Re      = 1;
 x0.dx      = 0.01;
 x0.dy      = 0.01;
 x0.dt      = 0.01;
 x0.time       = 10;
 
 % Set the intial guess == current vector
 x0.U_im_x  = myzero;
 x0.U_im_y  = myzero;
 % set solver parameters
 maxiters = 30;
 tol      = 1e-5;
 
% Preconditioning type
  pretype = 0;
  [x iter] = Newton_Jacobian_Appx('RHS_Newton_Test',pretype, x0, maxiters, tol); 
  %[x iter] = JFNK('RHS_Newton',pretype, x0, maxiters, tol); 
  
  U_im_x = x.U_im_x;
  U_im_y = x.U_im_y;
  
