%function [Ucat_phy_x, Ucat_phy_y, Pressure_phy, dx, dy, t] = Main(M, N)
function [Ucat_cal_x, Ucat_cal_y, Pressure_cal, dx, dy, t] = Main(M, N) % Debug

close all; format long;
ENABLE_VISUAL_GRID = 1;
ENABLE_CALCULATION = 0;
ENABLE_VISUAL_PLOT = 1;

%% Variables' dimension in non-staggered grid
%M = 8;
%N = 8;

%% Ghost variables' dimension in non-staggered grid
M2 = M+2;
N2 = N+2;

%% Variables' dimension in staggered grid
M3 = M+3;
N3 = N+3;

dt = 1E-4;
MAXTIME = 2;%round(1/dt);

%% Physical parameters
Re = 1;
L = 1;
U = 1;

%% Initialization process
dU_x = 0;
dU_y = 0;
t = 0;
[Ucont_x, Ucont_y, ...
    Ucat_phy_x, Ucat_phy_y, Pressure_phy, ...
    Ucat_cal_x, Ucat_cal_y, Pressure_cal, ...
    Ubcs_x, Ubcs_y, Pbcs, ...
    dx, dy] = Init(M, N, M2, N2, M3, N3, L, U, ENABLE_VISUAL_GRID)

enforce_bcs;

if ENABLE_VISUAL_GRID
    plot_coor_in_grid;
end

if ENABLE_VISUAL_PLOT
    myplot;
end

%% Calculation process
while ENABLE_CALCULATION
    for time_step = 1:MAXTIME       

        tic
        t = time_step * dt
    
        [U_im_x, U_im_y] =  Momentum_Solver(dU_x , dU_y, Ucont_x, Ucont_y, Ucat_cal_x, Ucat_cal_y, Ubcs_x, Ubcs_y, Pressure_cal, Re, dx, dy, dt, t);        
            
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