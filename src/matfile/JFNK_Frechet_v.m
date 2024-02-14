function F_v = JFNK_Frechet_v(F,u, v)

M = length(u.Ucont_x(:,1));
N = length(u.Ucont_x(1,:));
n= 2*M*N;

n = length(v);

% Current evalutation of function
P   = feval(F,u);
epsilon = 1e-9;

u_not = u;

[d_x d_y] = Un_Vectorize(v,M,N);


u_not.U_im_x = u.U_im_x + epsilon * d_x;
u_not.U_im_y = u.U_im_y + epsilon * d_y;

P_1 = feval(F,u_not);
F_v = (P_1 - P) / epsilon;


