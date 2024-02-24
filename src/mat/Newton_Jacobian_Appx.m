function [x its]= Newton_Jacobian_Appx(myfunc,pretype, x0, maxits, tol)
M = length(x0.Ucont_x(:,1));
N = length(x0.Ucont_x(1,:));
% Set the initial guess
initial_guess = [zeros(M*N,1); zeros(M*N,1)];

x = x0;
r0 = norm(feval(myfunc,x0),2);
e  = r0;

%%
% specify details for FGMRES solver
im = 30;        % restarted
max_subit = 60; % No more than 30 iterations  
tolIts  = tol; % only drops down by 5 order of magnitude
p  =  30;

for  iter = 1: maxits
  
  e   = feval(myfunc,x);
  Newton_Krylov_norm = norm(e,inf);
  
  fprintf('subiter=%d - NK norm = %f\n ',iter,Newton_Krylov_norm);
  if (Newton_Krylov_norm < tol*r0 || iter == maxits)
      its = iter;
      break;
  else
      

    %J_b = Momentum_Jacobian_Appx(myfunc,x);
    if (mod(iter,p) == 0 || iter == 1)
        J_b = Jacobian_Assemble(x);
        
            % Preconditioning techniques
         if pretype ==  0
         fprintf(' Newton with approximated J - No Precond \n');
         PRE = struct('L',sparse(0,0),'U',sparse(0,0));
         precfun = 'precNoPRE';
         end
        
        if pretype == 1
        
        %[L, U] = ilu0(A);
         fprintf(' Newton with approximated J - LU precond \n');
        ILU.L = spdiags(spdiags(J_b,0),0,2*M*N,2*M*N);%L;
        ILU.U = spdiags(spdiags(J_b,0),0,2*M*N,2*M*N);%U;
        PRE = ILU;
        precfun='precLU';
        end
   
        if pretype == 2
          fprintf(' Newton with approximated J - SOR precond \n');
        L = tril(J_b); D = diag(J_b); U = triu(J_b);
        PRE.L = L; 
        PRE.D = full(D);
        PRE.U = U;
        PRE.kiter = 1;
        precfun='precLDU';
        end
    
    end
    %figure(1)
    %spy(J_b,'b')
    %hold on
    %spy(J_c,'r')
    
    %x_b = diag(J_b,0);
    %x_c = diag(J_c,0);
    %figure(3)
    %plot(x_b,'b x');
    %hold on
    %plot(x_c,'r o');
    %pause
    %norm(J_c-J_b,inf)
       
    A   = J_b;
    rhs = -e; 
    
    
    % Solve to advance the solution
    d_x = fgmres (A,PRE,precfun,rhs,initial_guess,im,max_subit,tolIts);
    
    [delta_U delta_V] = Un_Vectorize(d_x,M,N);
    
    % Remember
    if (iter == 1)
        dU0 = delta_U;
        dV0 = delta_V;
    end
    % Advance the solution
    % apply line search technique    
    neta0 = 0.2;
    neta1 = 0.5;
    
    rU = (neta1*delta_U + neta0*dU0);
    rV = (neta1*delta_V + neta0*dV0);    
    
    x.U_im_x   = x.U_im_x + rU;
    x.U_im_y   = x.U_im_y + rV;
    
    dU0 = rU;
    dV0 = rV;
        
  end
end