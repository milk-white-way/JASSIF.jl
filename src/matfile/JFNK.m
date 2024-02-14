function [x its]= JFNK(myfunc,pretype,PRE,precfun, x0, maxits, tol)
M = length(x0.Ucont_x(:,1));
N = length(x0.Ucont_x(1,:));
zero_guess = [zeros(M*N,1); zeros(M*N,1)];
x = x0;
r0 = norm(feval(myfunc,x0),2);
e  = r0;

% specify details for FGMRES solver
im = 60;        % restarted
max_subit = 60; % No more than 30 iterations  
tolIts  = tol; % only drops down by 5 order of magnitude
p = maxits;         % Evaluate the Jacobian after p Newton iteration

for  iter = 1: maxits
  
  e   = feval(myfunc,x);
  
  if (norm(e,inf) < tol*r0 || iter ==maxits)
      its = iter;
      break;
  else
    
       
    fprintf('\nsubiter= %d - NK norm = %f \n ',iter,norm(e,inf));
    rhs = -e;
    
    %Preconditioning technique---------------    
    if (mod(iter,p) == 0||iter==1)
    
    if pretype == 0
          fprintf('JFNK with no precond \n');
    end
        
    if pretype == 1
     
     fprintf('JFNK with ILU precond \n');
     A = Jacobian_Assemble(x0); 
      
     OPTS.droptol = 0.05;
     OPTS.thresh  = 0;      
     %[L, U] = luinc(A,OPTS);%ilu0(A);
     ILU.L = spdiags(sqrt(spdiags(A,0)),0,2*M*N);%L;
     ILU.U = spdiags(sqrt(spdiags(A,0)),0,2*M*N);%U;
     PRE = ILU;
     precfun='precLU';
    end
   
    if pretype == 2
     
     A = Jacobian_Assemble(x0); 
     fprintf('JFNK with SOR precond \n');
     L = tril(A); D = diag(A); U = triu(A);
     PRE.L = L; 
     PRE.D = full(D);
     PRE.U = U;
     PRE.kiter = 2;
     precfun='precLDU';
    end
    end
    % End of preconditioning -----------------
    
    % Solve to advance the solution
    d_x = fgmres_JFNK_Frechet (myfunc,PRE,precfun,rhs,x,im,max_subit,tolIts);
    
    [delta_U delta_V] = Un_Vectorize(d_x,M,N);
    
    % Advance the solution
    x.U_im_x   = x.U_im_x + delta_U;
    x.U_im_y   = x.U_im_y + delta_V;
    
    
  end
end