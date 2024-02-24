function [RHS_vector] = RHS_Newton(OPTS)


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
     
          % Forming BCSs & convert to cartesian components
          [Ucat_x Ucat_y U_im_x U_im_y] = FormBCS(U_im_x, U_im_y, Ubcs_x, Ubcs_y,dx,dy,Re,time);    
 
          
          [RHS_x RHS_y] = RHS_Calculation(U_im_x, U_im_y, Ucat_x, Ucat_y, Pressure,Re,dx,dy);

          RHS_x = RHS_x -(1.5/dt) * (U_im_x  - Ucont_x) + (0.5/dt) * dU_old_x;
          RHS_y = RHS_y -(1.5/dt) * (U_im_y  - Ucont_y) + (0.5/dt) * dU_old_y;  
          
          M = length(Ucont_x(:,1));
          N = length(Ucont_x(1,:));

          for j= 1:N
              RHS_x(1,j) = 0;
              RHS_x(M-1,j) = 0;
              RHS_x(M,j)  =0;
              
          end
          
          for i =1:M
              RHS_y(i,1) = 0;
              RHS_y(i,N) = 0;
              RHS_y(i,N-1) = 0;
          end
          
          v_1 = Vectorize(RHS_x);
          v_2 = Vectorize(RHS_y);
          
          RHS_vector = [v_1; v_2];
         
 