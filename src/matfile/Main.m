function [Ucat_x Ucat_y Pressure dx dy M N t] = Main()

clear all
format long;

M = 3;
N = 3;
% Add two more ghost points
M2 = M+2;
N2 = N+2;

dt = 0.05;
MAXTIME = 200;%round(1/dt);
%Note for Taylor Problem Re=1
Re = 100;
% Remember dU for RK
dU_x = 0;
dU_y = 0;
% Initialization process
[Ucont_x Ucont_y Ucat_x Ucat_y Ubcs_x Ubcs_y Pressure dx dy] = Init(M2,N2);%_Taylor_Green(M2,N2);

% It is the time integration scheme
for time_step = 1:MAXTIME       
       
    tic
    time_step
    t = time_step * dt    
       
    [U_im_x U_im_y] =  Momentum_Solver(dU_x , dU_y, Ucont_x, Ucont_y, Ucat_x, Ucat_y, Ubcs_x, Ubcs_y, Pressure, Re, dx, dy, dt,t);        
     
    [phi]  = Poisson_Solver(U_im_x, U_im_y, dx, dy, dt);
    
    [U_x_new U_y_new P_new] = Update_Solution(U_im_x, U_im_y, Pressure, phi,dx,dy,dt);
    
    %Update dU
    dU_x = U_x_new - Ucont_x;
    dU_y = U_y_new - Ucont_y;
    
    % Go next time step
    [Ucat_x Ucat_y Ucont_x Ucont_y] = FormBCS(U_x_new, U_y_new, Ubcs_x, Ubcs_y, dx, dy,Re,t);
    
    Pressure = P_new;          
        
    % Assess the divergence of flow field
    Div = Divergence(U_x_new, U_y_new,dx ,dy);
    MaxDiv = norm(Div,inf)
   
    toc             
    
    norm(Vectorize(dU_x),inf)
end
    

