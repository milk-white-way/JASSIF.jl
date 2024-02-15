function [Ucat_phy_x, Ucat_phy_y, Pressure_phy, dx, dy, t] = Main()

clearvars; close all;

format long;

% Variables' dimension in non-staggered grid
M = 4;
N = 4;

% Ghost variables' dimension in non-staggered grid
M2 = M+2;
N2 = N+2;

% Variables' dimension in staggered grid
M3 = M+3;
N3 = N+3;

dt = 0.05;
MAXTIME = 200;%round(1/dt);

Re = 100;

% Initialization process
dU_x = 0;
dU_y = 0;
t = 0;
[Ucont_x, Ucont_y, ...
    Ucat_phy_x, Ucat_phy_y, Pressure_phy, ...
    Ucat_cal_x, Ucat_cal_y, Pressure_cal, ...
    Ubcs_x, Ubcs_y, Pbcs, dx, dy] = Init(M,N, M2, N2, M3, N3)

% It is the time integration scheme
CALCULATE = 0;
while CALCULATE
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
    if time_step == MAXTIME
        CALCULATE = 0;
    end
end
    

