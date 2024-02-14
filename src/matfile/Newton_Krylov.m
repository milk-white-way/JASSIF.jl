function [x its]= Newton_Krylov(myfunc,pretype, x0, maxits, tol)

x = x0;
r0 = norm(feval(myfunc,x0),2);
e  = r0;
n = length(x0);
% specify details for FGMRES solver
im = 30;        % restarted
max_subit = 60; % No more than 30 iterations  
tolIts  = tol; % only drops down by 5 order of magnitude


for  iter = 1: maxits
  
  e   = feval(myfunc,x);
  
  if norm(e,inf) < tol*r0
      its = iter;
      break;
  else
    J_b = Jacobian_Appx(myfunc,x);
    A   = J_b;
    rhs = -e; 
    
    % Preconditioning techniques
     if pretype ==  0
      PRE = struct('L',sparse(0,0),'U',sparse(0,0));
      precfun = 'precNoPRE';
     end
    
     if pretype == 1
        
     [L, U] = ilu0(A);
     ILU.L = L;
     ILU.U = U;
     PRE = ILU;
     precfun='precLU';
     end
   
    if pretype == 2
     L = tril(A); D = diag(A); U = triu(A);
     PRE.L = L; 
     PRE.D = full(D);
     PRE.U = U;
     PRE.kiter = 1;
      precfun='precLDU';
    end
    
    % Solve to advance the solution
    d_x = fgmres (A,PRE,precfun,rhs,zeros(n,1),im,max_subit,tolIts)
    x   = x + d_x;
    
  end
end