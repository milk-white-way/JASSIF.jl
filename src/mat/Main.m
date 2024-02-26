function [PhysDom, CompDom, HaloDom, FluxSum, dx, dy, t] = ...
    Main(M, N) % Get rid of HaloDom for production run

    close all; 
    format shortE;
    format compact;

    IMPORT_HDF5 = 'PostMomentumData.h5';
    if exist(IMPORT_HDF5, 'file') == 2
        delete(IMPORT_HDF5);
    end

    AMRESSIF = 'input';
    if exist(AMRESSIF, 'file') == 2
        delete(AMRESSIF);
    end

    %% Some flags
    ENABLE_VISUAL_GRID = 0;
    ENABLE_CALCULATION = 1;
    ENABLE_VISUAL_PLOT = 0;
    ENABLE_BC_PERIODIC = 1;

    ENABLE_DEBUGGING = 0;

    Nghost = 2;

    %% Physical parameters
    Re = 1;
    L = 1;
    U = 1;
    %% Solver parameters
    dt = 1E-4;
    MAXTIME = 1;
    %% AMRESSIF parameters
    max_grid_size = 8;
    plot_int = 1;

    if ENABLE_BC_PERIODIC
        bc_lo = [0, 0];
        bc_hi = [0, 0];
    else
        error('Non-periodic boundary condition is not supported yet!')
    end

    fprintf('\tJust Another Simulation Suite For Incompressible Flows \n');
    fprintf('\t \t \t \t \tVersion 3.0 \n');
    fprintf('\t \t \tAuthors: Trung Le and Tam Nguyen \n');
    fprintf('\t \t2022-2024 North Dakota State University \n');
    fprintf('\n=================================================================\n');

    set(0,'DefaultFigureWindowStyle','docked')
    %% Variables' dimension in non-staggered grid
    %M = 8; % Now taking input
    %N = 8; % Now taking input
    fprintf('INFO: \tNumber of cells in x-direction = %d \n', M);
    fprintf('INFO: \tNumber of cells in y-direction = %d \n', N);
    fprintf('INFO: \tNumber of ghost layer = %d \n', Nghost);
    fprintf('INFO: \tTime step = %f \n', dt);

    %% Initialization process
    dU_x = 0;
    dU_y = 0;
    t = 0;

    tic;
    fprintf('INFO: \tBegin Initialization... ');
    [PhysDom, CompDom, FluxSum, M2, N2, M3, N3, iphys, iphye, jphys, jphye, dx, dy] = Init(M, N, Nghost, L, U, ENABLE_VISUAL_GRID);
    fprintf('Done! \n');
    fprintf('INFO: \tSpatial discretization: dx = %f, dy = %f \n', dx, dy);

    [CompDom, HaloDom, ~, ~, ~, ~, ~] = TAM_enforce_bcs_v2(PhysDom, CompDom, M, N, M2, N2, Nghost, iphys, iphye, jphys, jphye, ENABLE_BC_PERIODIC, ENABLE_DEBUGGING);

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

    % Optional step: Create input file for AMRESSIF Poisson solver
    fid = fopen(AMRESSIF, 'w');
        fprintf(fid, '#Hello, from JASSIF! \n');
        fprintf(fid, 'n_cell = %d \n', M);
        fprintf(fid, 'max_grid_size = %d \n', max_grid_size);
        fprintf(fid, 'plot_int = %d \n', plot_int);
        fprintf(fid, 'time = %d \n', t);
        fprintf(fid, 'bc_lo = %d %d \n', bc_lo(1), bc_lo(2));
        fprintf(fid, 'bc_hi = %d %d \n', bc_hi(1), bc_hi(2));

    % Close the file
    fclose(fid);

    %% Calculation process
    while ENABLE_CALCULATION
        fprintf('\nINFO: \tBegin Calculation... \n');
        for time_step = 1:MAXTIME       

            tic;
            t = time_step * dt;

            %% Solve the divergence-free Momentum Equation to obtain next-timestep contravariant velocity components
            [FluxSum, Ucont_im_x, Ucont_im_y] = Runge_Kutta(CompDom, FluxSum, dU_x, dU_y, M, N, M2, N2, M3, N3, Nghost, iphys, iphye, jphys, jphye, Re, dx, dy, dt, t, ENABLE_BC_PERIODIC, ENABLE_DEBUGGING);

            %{
            %% Solve the Poisson Equation to obtain correction field 'phi'
            % Step 1: Export contravariant velocity components and pressure field to hdf5 file
            h5create(IMPORT_HDF5, '/Ucont/imx', [N M3]);
            h5create(IMPORT_HDF5, '/Ucont/imy', [N3 M]);
            h5create(IMPORT_HDF5, '/UserCtx/p', [N M]);

            h5write(IMPORT_HDF5, '/Ucont/imx', Ucont_im_x);
            h5write(IMPORT_HDF5, '/Ucont/imy', Ucont_im_y);
            h5write(IMPORT_HDF5, '/UserCtx/p', PhysDom.Pressure);

            % Step 2: Update time in AMRESSIF input file
            fid = fopen(AMRESSIF);
            S = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
            S = S{1};
            idx = contains(S, 'time');
            S{idx} = ['time = ' num2str(t)];
            fid = fopen(AMRESSIF, 'w');
            fprintf(fid, '%s\n', S{:});
            fclose(fid);

            % Step 2: Call in the Poisson solver from AMRESSIF source code to solve the Poisson Equation

            %}

            % Calculate the divergence of the new contraction velocity field
            P_Div_Phys = Divergence(Ucont_im_x, Ucont_im_y, M, N, dx, dy);
            % Ensemble the right-hand side of the Poisson equation with one layer of ghost cells
            P_Div = zeros(N2, M2);

            b = zeros(M*N, 1);
            for ii = 1:M
                for jj = 1:N
                    index = glidx(ii, jj, M);

                    b(index) = P_Div(jj, ii);
                end
            end

            A = Poisson_LHS_Neumann(M, N, dx, dy);
            phi = A\b';

            %[phi] = Poisson_Solver(U_im_x, U_im_y, dx, dy, dt);

            %{
            [U_x_new, U_y_new, P_new] = Update_Solution(U_im_x, U_im_y, Pressure, phi, dx, dy, dt);
                
            %Update dU
            dU_x = U_x_new - Ucont_x;
            dU_y = U_y_new - Ucont_y;
                
            % Go next time step
            [Ucat_x, Ucat_y, Ucont_x, Ucont_y] = FormBCS(U_x_new, U_y_new, Ubcs_x, Ubcs_y, dx, dy,Re,t);
                
            Pressure = P_new;          
                
            % Assess the divergence of flow field
            Div = Divergence(U_x_new, U_y_new,dx ,dy);
            MaxDiv = norm(Div, inf)
            
            norm(Vectorize(dU_x), inf)
            %}
            fprintf('INFO: \t Time step %d is done! \n', time_step);
            toc;
        end

        if time_step == MAXTIME
            ENABLE_CALCULATION = 0;
        end
    end
end