function [RHS_vector] = RHS_Newton_Test(OPTS)


dU_old_x = OPTS.dU_x;
dU_old_y = OPTS.dU_y;
Ucont_x  = OPTS.Ucont_x;
Ucont_y  = OPTS.Ucont_y;
U_im_x   = OPTS.U_im_x;
U_im_y   = OPTS.U_im_y;

Ucat_x   = OPTS.Ucat_x;
Ucat_y   = OPTS.Ucat_y;
Ubcs_x   = OPTS.Ubcs_x;
Ubcs_y   = OPTS.Ubcs_y;
Pressure = OPTS.Pressure;
Re       = OPTS.Re;
dx       = OPTS.dx;
dy       = OPTS.dy;
dt       = OPTS.dt;
time     = OPTS.time;
     
M = length(Ucont_x(:,1));
N = length(Ucont_x(1,:));
  

for i = 1:M
    for j=1:N
   RHS_x(i,j) =U_im_x(i,j).^3- 2;
    end
end

for i = 1:M
    for j=1:N
   RHS_y(i,j) = U_im_y(i,j).^3 -2;
    end
end

          v_1 = Vectorize(RHS_x);
          v_2 = Vectorize(RHS_y);
          
          
          RHS_vector = [v_1; v_2];
         
 