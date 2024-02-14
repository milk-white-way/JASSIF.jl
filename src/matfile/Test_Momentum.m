function Main()

clear all
MAXTIME =1000;
M = 80;
N = 80;
% Add two more ghost points
M2 = M+2;
N2 = N+2;

dt = 1e-3;
Re = 1;
% Remember dU for RK
dU_x = 0;
dU_y = 0;
% Initialization process
[Ucont_x Ucont_y Ucat_x Ucat_y Ubcs_x Ubcs_y Pressure dx dy] = Init_Taylor_Green(M2,N2);

t = 0;
    
[U_im_x U_im_y] =  Runge_Kutta(dU_x , dU_y, Ucont_x, Ucont_y, Ucat_x, Ucat_y, Ubcs_x, Ubcs_y, Pressure, Re, dx, dy, dt,t);            
     
%[RHS_x RHS_y] =  Momentum_RHS(dU_x , dU_y, Ucont_x,Ucont_y,Ucont_x, Ucont_y, Ucat_x, Ucat_y, Ubcs_x, Ubcs_y, Pressure, Re, dx, dy, dt);            

   
figure(1)
%mesh(U_im_x- Ucont_x)
mesh(Ucat_x)
