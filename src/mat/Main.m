function Main(control_path, checkpoint_path)

    format shortE;
    format compact;

    %% Include Control Parameters
    run(control_path);
    checkpoint_form = num2str(floor(log10(MAXTIME)) + 1);
    dir_name = '/Checkpoints/';
    full_path = strcat(checkpoint_path, dir_name);
    if ~exist(full_path, 'dir')
        mkdir(full_path);
    end
    chkname = strcat('%schk_%0', checkpoint_form, 'd.mat');

    %% Extra solver's parameters that for dev
    Nghost = 2;

    %% AMRESSIF parameters
    if ENABLE_AMRESSIF
        IMPORT_HDF5 = 'PostMomentumData.h5';
        if exist(IMPORT_HDF5, 'file') == 2
            delete(IMPORT_HDF5);
        end

        AMRESSIF = 'input';
        if exist(AMRESSIF, 'file') == 2
            delete(AMRESSIF);
        end

        max_grid_size = 8;
        plot_int = 1;

        if ENABLE_BC_PERIODIC
            bc_lo = [0, 0];
            bc_hi = [0, 0];
        else
            error('Non-periodic boundary condition is not supported yet!')
        end
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

    %% Writing initial configuration as checkpoints
    full_name = sprintf(chkname, full_path, t);
    save(full_name, 'PhysDom', 'FluxSum', 'dx', 'dy', 't')

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

    if ENABLE_AMRESSIF
        %% Write the input file for AMRESSIF
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
    end

    %% Calculation process
    while ENABLE_CALCULATION
        tic;
        fprintf('\nINFO: \tBegin Calculation... \n');
        for time_step = 1:MAXTIME       

            t = time_step * dt;

            Ucont_pre_x = CompDom.Ucont_x;
            Ucont_pre_y = CompDom.Ucont_y;

            %% Solve the divergence-free Momentum Equation to obtain next-timestep contravariant velocity components
            [FluxSum, Ucont_im_x, Ucont_im_y] = Runge_Kutta(CompDom, FluxSum, dU_x, dU_y, M, N, M2, N2, M3, N3, Nghost, iphys, iphye, jphys, jphye, Re, dx, dy, dt, t, ENABLE_BC_PERIODIC, ENABLE_DEBUGGING);

            if ENABLE_AMRESSIF
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
            end

            % Calculate the divergence of the new contraction velocity field
            P_Div_Phys = Divergence(Ucont_im_x, Ucont_im_y, M, N, dx, dy);
            % Ensemble the right-hand side of the Poisson equation with one layer of ghost cells
            M4 = M3 + 1;
            N4 = N3 + 1;
            P_Div = zeros(N4, M4);
            P_Div(2:N3, 2:M3) = P_Div_Phys;

            % Apply periodic boundary condition
            P_Div(2:N3, 1) = P_Div_Phys(:, end);
            P_Div(2:N3, end) = P_Div_Phys(:, 1);
            P_Div(1, 2:M3) = P_Div_Phys(end, :);
            P_Div(end, 2:M3) = P_Div_Phys(1, :);

            % Ensemble the RHS vector
            b = zeros(N4*M4, 1);
            for ii = 1:M4
                for jj = 1:N4
                    index = glidx(ii, jj, M4, N4);

                    b(index) = P_Div(jj, ii);
                end
            end

            % Ensemble the LHS sparse matrix
            A = Poisson_LHS_Neumann(M4, N4, dx, dy);
            A = A';

            % Solve the Poisson equation
            phi = A\b;

            % Return the solution to the physical domain
            for index = 1:M4*N4
                [ii, jj] = lidx(index, M4);
                a_Phi(jj, ii) = phi(index);
            end

            % Calculate the gradient of phi
            for ii = 1:M3
                for jj = 1:N
                    PGrad_x(jj, ii) = (a_Phi(jj, ii+1) - a_Phi(jj, ii))/dx;
                end
            end

            for ii = 1:M
                for jj = 1:N3
                    PGrad_y(jj, ii) = (a_Phi(jj+1, ii) - a_Phi(jj, ii))/dy;
                end
            end

            % Update corrected contravariant velocity components
            Ucont_new_x = Ucont_im_x - PGrad_x*dt/1.5;
            Ucont_new_y = Ucont_im_y - PGrad_y*dt/1.5;

            % Update dU
            dU_x = Ucont_new_x - Ucont_pre_x;
            dU_y = Ucont_new_y - Ucont_pre_y;

            %% Update all new fields to physical domain
            % Update pressure field
            for ii = 1:M
                for jj = 1:N
                    PhysDom.Pressure(jj, ii) = PhysDom.Pressure(jj, ii) + a_Phi(jj+1, ii+1);
                end
            end

            % Update Carterian velocity components
            [Ucat_new_x, Ucat_new_y] = Contra_To_Cart(Ucont_new_x, Ucont_new_y, M, N);
            PhysDom.Ucat_x = Ucat_new_x;
            PhysDom.Ucat_y = Ucat_new_y;

            %% Ensemble the computational domain before next time step
            CompDom.Ucont_x = Ucont_new_x;
            CompDom.Ucont_y = Ucont_new_y;

            %% Enforce boundary conditions
            [CompDom, ~, ~, ~, ~, ~, ~] = TAM_enforce_bcs_v2(PhysDom, CompDom, M, N, M2, N2, Nghost, iphys, iphye, jphys, jphye, ENABLE_BC_PERIODIC, ENABLE_DEBUGGING);

            % Assess the divergence of flow field
            Div = Divergence(Ucont_new_x, Ucont_new_y, M, N, dx , dy);
            MaxDiv = norm(Div, inf);
            
            norm(Vectorize(dU_x), inf);
            fprintf('INFO: \t Time Step %d is done where the time is %.4f \n', time_step, t);

            %% Writing checkpoints
            if rem(time_step, checkpoint_freq) == 0
                full_name = sprintf(chkname, full_path, time_step);
                save(full_name, 'PhysDom', 'FluxSum', 'dx', 'dy', 't')
            end
        end

        %% Check the time
        if time_step == MAXTIME
            ENABLE_CALCULATION = 0;
        end
    end
    toc;

    %% Post-processing: plot the final results
    if ENABLE_VISUAL_PLOT
        TAM_myplot();
    end
end