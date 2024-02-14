function [U_im_x U_im_y] =  Momentum_Solver(dU_x , dU_y, Ucont_x, Ucont_y, Ucat_x, Ucat_y, Ubcs_x, Ubcs_y, Pressure, Re, dx, dy, dt,t);        
 
 % Set initial guess
 x0.dU_x    = dU_x;
 x0.dU_y    = dU_y;
 x0.Ucont_x = Ucont_x;
 x0.Ucont_y = Ucont_y;
 x0.Ucat_x  = Ucat_x;
 x0.Ucat_y  = Ucat_y;
 x0.Ubcs_x  = Ubcs_x;
 x0.Ubcs_y  = Ubcs_y;
 x0.Pressure= Pressure;
 x0.Re      = Re;
 x0.dx      = dx;
 x0.dy      = dy;
 x0.dt      = dt;
 x0.time       = t;
 M = length(Ucont_x(:,1));
 N = length(Ucont_x(1,:));
 Momentum = 0;
 
 
 % Set the intial guess == current vector
 x0.U_im_x  = Ucont_x;
 x0.U_im_y  = Ucont_y;
 % set solver parameters
 maxiters = 100;
 tol      = 1e-5;
 pretype = 0;


fprintf('Momentum solver \n');
 % Preconditioning techniques
    if pretype ==  0
       A = speye(2*M*N);
       PRE = struct('L',sparse(0,0),'U',sparse(0,0));
       precfun = 'precNoPre';
    end
    
    if pretype == 1
        
     A = Jacobian_Assemble(x0); 
      
     OPTS.droptol = 0.05;
     OPTS.thresh  = 0;      
     %[L, U] = luinc(A,OPTS);%ilu0(A);
     ILU.L = spdiags(sqrt(spdiags(A,0)),0,2*M*N,2*M*N);%L;
     ILU.U = spdiags(sqrt(spdiags(A,0)),0,2*M*N,2*M*N);%U;
     PRE = ILU;
     precfun='precLU';
    end
   
    if pretype == 2
     
     A = Jacobian_Assemble(x0); 
     
     L = tril(A); D = diag(A); U = triu(A);
     PRE.L = L; 
     PRE.D = full(D);
     PRE.U = U;
     PRE.kiter = 2;
     precfun='precLDU';
    end
    
    
  fprintf('Complete Analytical Jacobian...\n');  
    
 %------------------  Momentum solvers -------------------------
 if Momentum == 0
   [U_im_x U_im_y] =  Runge_Kutta(dU_x , dU_y, Ucont_x, Ucont_y, Ucat_x, Ucat_y, Ubcs_x, Ubcs_y, Pressure, Re, dx, dy, dt,t);        
 end
 
 if Momentum == 1
  [x iter] = Newton_Jacobian_Appx('RHS_Newton',pretype, x0, maxiters, tol); 
  U_im_x = x.U_im_x;
  U_im_y = x.U_im_y;
  
 end
 
 if Momentum == 2
  [x iter] = JFNK('RHS_Newton',pretype,PRE,precfun, x0, maxiters, tol); 
  U_im_x = x.U_im_x;
  U_im_y = x.U_im_y;
  
 end

 