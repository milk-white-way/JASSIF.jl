function [PhysDom, CompDom, HaloDom, dx, dy, t] = Main(M, N) % Get rid of HaloDom for production run
    close all; format long;
    %% Some flags
    ENABLE_VISUAL_GRID = 1;
    ENABLE_CALCULATION = 0;
    ENABLE_VISUAL_PLOT = 0;
    ENABLE_BC_PERIODIC = 1;

    ENABLE_DEBUGGING = 1;

    %% Physical parameters
    Re = 1;
    L = 1;
    U = 1;
    %% Solver parameters
    dt = 1E-4;
    MAXTIME = 2;

    fprintf('\tJust Another Simulation Suite For Incompressible Flows \n');
    fprintf('\t \t \t \t \tVersion M-2024.2 \n');
    fprintf('\t \t \tAuthors: Trung Le and Tam Nguyen \n');
    fprintf('\n=================================================================\n');

    set(0,'DefaultFigureWindowStyle','docked')
    %% Variables' dimension in non-staggered grid
    %M = 8;
    %N = 8;
    fprintf('INFO: \tNumber of cells in x-direction = %d \n', M);
    fprintf('INFO: \tNumber of cells in y-direction = %d \n', N);

    Nghost = 2; fprintf('INFO: \tNumber of ghost layer = %d \n', Nghost);

    %% Initialization process
    dU_x = 0;
    dU_y = 0;
    t = 0;
    fprintf('INFO: \tTime step = %f \n', dt);

    tic;
    fprintf('INFO: \tBegin Initialization... ');
    [PhysDom, CompDom, HaloDom, M2, N2, dx, dy] = Init(M, N, Nghost, L, U, ENABLE_VISUAL_GRID);
    fprintf('Done! \n');
    fprintf('INFO: \tSpatial discretization: dx = %f, dy = %f \n', dx, dy);

    TAM_enforce_bcs();

    if ENABLE_VISUAL_GRID
        fprintf('INFO: \tVisualizing grid... ');
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
        fprintf('\nINFO: \tBegin Calculation... \n');
        for time_step = 1:MAXTIME       

            tic;
            t = time_step * dt;
        
            [FluxFld, U_im_x, U_im_y] =  Runge_Kutta(dU_x, dU_y, Ucont_x, Ucont_y, Ucat_cal_x, Ucat_cal_y, Ubcs_x, Ubcs_y, Pressure_cal, Re, dx, dy, dt, t);        
                
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