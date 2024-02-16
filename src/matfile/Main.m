%function [Ucat_phy_x, Ucat_phy_y, Pressure_phy, dx, dy, t] = Main(M, N)
function [Ucat_cal_x, Ucat_cal_y, Pressure_cal, dx, dy, t] = Main(M, N) % Debug
    close all; format long;
    %% Some flags
    ENABLE_VISUAL_GRID = 0;
    ENABLE_CALCULATION = 0;
    ENABLE_VISUAL_PLOT = 1;
    ENABLE_BC_PERIODIC = 1;

    ENABLE_DEBUGGING = 0;

    %% Physical parameters
    Re = 1;
    L = 1;
    U = 1;
    %% Solver parameters
    dt = 1E-4;
    MAXTIME = 2;

    fprintf('\t Just Another Simulation Suite For Incompressible Flows \n');
    fprintf('\t \t \t \t \t Version M-2024.2 \n');
    fprintf('\t \t \t Authors: Trung Le and Tam Nguyen \n');
    fprintf('\n=================================================================\n');

    %% Variables' dimension in non-staggered grid
    %M = 8;
    %N = 8;
    fprintf('INFO: \t Number of cells in x-direction = %d \n', M);
    fprintf('INFO: \t Number of cells in y-direction = %d \n', N);

    %% Ghost variables' dimension in non-staggered grid
    Nghost = 1; fprintf('INFO: \t Number of ghost layer = %d \n', Nghost);
    M2 = M + (2*Nghost);
    N2 = N + (2*Nghost);

    %% Variables' dimension in staggered grid
    M3 = M + (2*Nghost) + 1;
    N3 = N + (2*Nghost) + 1;

    %% Initialization process
    dU_x = 0;
    dU_y = 0;
    t = 0;
    fprintf('INFO: \t Time step = %f \n', dt);

    tic;
    fprintf('INFO: \t Begin Initialization... ');
    [Ucont_x, Ucont_y, ...
        Ucat_phy_x, Ucat_phy_y, Pressure_phy, ...
        Ucat_cal_x, Ucat_cal_y, Pressure_cal, ...
        Ubcs_x, Ubcs_y, Pbcs, ...
        dx, dy] = Init(M, N, M2, N2, M3, N3, L, U, ENABLE_VISUAL_GRID);
    fprintf('Done! \n');

    TAM_enforce_bcs();

    if ENABLE_VISUAL_GRID
        fprintf('INFO: \t Visualizing grid... ');
        TAM_coordinates();
        H = gcf;
        fprintf('Done! Computational grid is represented in Figure %d. \n', H.Number);
    end

    if ENABLE_VISUAL_PLOT
        TAM_myplot();
    end
    toc;

    %% Calculation process
    while ENABLE_CALCULATION
        fprintf('\nINFO: \t Begin Calculation... \n');
        for time_step = 1:MAXTIME       

            tic;
            t = time_step * dt;
        
            [U_im_x, U_im_y] =  Runge_Kutta(dU_x, dU_y, Ucont_x, Ucont_y, Ucat_cal_x, Ucat_cal_y, Ubcs_x, Ubcs_y, Pressure_cal, Re, dx, dy, dt, t);        
                
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
            MaxDiv = norm(Div, inf)
            
            norm(Vectorize(dU_x), inf)

            fprintf('INFO: \t Time step %d is done! \n', time_step);
            toc;
        end

        if time_step == MAXTIME
            CALCULATE = 0;
        end
    end
end