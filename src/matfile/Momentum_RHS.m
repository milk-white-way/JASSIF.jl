function [RHS_x RHS_y] = Momentum_RHS(dU_old_x, dU_old_y, U_im_x, U_im_y, Ucont_x, Ucont_y, Ucat_x, Ucat_y, Ubcs_x, Ubcs_y, Pressure, Re,dx,dy, dt)

[RHS_x RHS_y] = RHS_Calculation(U_im_x, U_im_y, Ucat_x, Ucat_y, Pressure,Re,dx,dy);

RHS_x = RHS_x -(1.5/dt) * (U_im_x  - Ucont_x) + (0.5/dt) * dU_old_x;
RHS_y = RHS_y -(1.5/dt) * (U_im_y  - Ucont_y) + (0.5/dt) * dU_old_y;     
      
