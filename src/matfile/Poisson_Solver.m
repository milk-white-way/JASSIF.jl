function [phi] = Poisson_Solver(Ucont_x, Ucont_y, dx, dy, dt)

M = length(Ucont_x(:,1));
N = length(Ucont_x(1,:));

% Assemble Matrix A
A = Poisson_LHS_Neumann(M,N,dx,dy);

% Assemble rhs - b
b = Poisson_RHS_Neumann(Ucont_x, Ucont_y, dx, dy, dt);

%phi = gmres(A,b',60,1e-8);

 phi = A\b';
