fprintf('\nINFO: \tEnforcing boundary conditions... \n'); 

if ENABLE_BC_PERIODIC
    %% North and South boundaries:
    % North wall is the same as the South wall
    PhysDom.Ubcs_x(1+Nghost:N+Nghost, Nghost+1) = PhysDom.Ucat_x(:, M);
    PhysDom.Ubcs_y(1+Nghost:N+Nghost, Nghost+1) = PhysDom.Ucat_y(:, M);
    PhysDom.Ubcs_x(1+Nghost:N+Nghost, M+Nghost) = PhysDom.Ucat_x(:, 1);
    PhysDom.Ubcs_y(1+Nghost:N+Nghost, M+Nghost) = PhysDom.Ucat_y(:, 1);

    % Halo region
    for ii = ( 1+Nghost ):( M+Nghost )
        for jj = 1:Nghost 
            mirror_dis = 1 + Nghost - jj;
            mirror_loc = N + Nghost - mirror_dis;
            %mirror_loc = Nghost + 1 + mirror_dis;

            HaloDom.Ucat_x(jj, ii) = CompDom.Ucat_x(mirror_loc, ii);
            HaloDom.Ucat_y(jj, ii) = CompDom.Ucat_y(mirror_loc, ii);
            HaloDom.Pressure(jj, ii) = CompDom.Pressure(mirror_loc, ii);

            if ENABLE_DEBUGGING
                fprintf('DEBUG: \tGhost node (%d, %d) ', ii, jj);
                fprintf('mirrors inner node (%d, %d)\n', ii, mirror_loc);
            end
        end

        for jj = (N+Nghost+1):N2
            mirror_dis = jj - (N + Nghost);
            mirror_loc = 1 + Nghost + mirror_dis;
            %mirror_loc = N + Nghost + mirror_dis;

            HaloDom.Ucat_x(jj, ii) = CompDom.Ucat_x(mirror_loc, ii);
            HaloDom.Ucat_y(jj, ii) = CompDom.Ucat_y(mirror_loc, ii);
            HaloDom.Pressure(jj, ii) = CompDom.Pressure(mirror_loc, ii);

            if ENABLE_DEBUGGING
                fprintf('DEBUG: \tGhost node (%d, %d) ', ii, jj);
                fprintf('mirrors inner node (%d, %d)\n', ii, mirror_loc);
            end
        end
    end

    %% East and West boundaries:
    % East wall is the same as the West wall
    PhysDom.Ubcs_x(Nghost+1, 1+Nghost:M+Nghost) = PhysDom.Ucat_x(N, :);
    PhysDom.Ubcs_y(Nghost+1, 1+Nghost:M+Nghost) = PhysDom.Ucat_y(N, :);
    PhysDom.Ubcs_x(N+Nghost, 1+Nghost:M+Nghost) = PhysDom.Ucat_x(1, :);
    PhysDom.Ubcs_y(N+Nghost, 1+Nghost:M+Nghost) = PhysDom.Ucat_y(1, :);

    for jj = ( 1+Nghost ):( N+Nghost )
        for ii = 1:Nghost
            mirror_dis = Nghost + 1 - ii;
            mirror_loc = M + Nghost - mirror_dis;

            HaloDom.Ucat_x(jj, ii) = CompDom.Ucat_x(jj, mirror_loc);
            HaloDom.Ucat_y(jj, ii) = CompDom.Ucat_y(jj, mirror_loc);
            HaloDom.Pressure(jj, ii) = CompDom.Pressure(jj, mirror_loc);

            if ENABLE_DEBUGGING
                fprintf('DEBUG: \tGhost node (%d, %d) ', ii, jj);
                fprintf('mirrors inner node (%d, %d)\n', mirror_loc, jj);
            end
        end

        for ii = (M+Nghost+1):M2
            mirror_dis = ii - (M + Nghost);
            mirror_loc = 1 + Nghost + mirror_dis;

            HaloDom.Ucat_x(jj, ii) = CompDom.Ucat_x(jj, mirror_loc);
            HaloDom.Ucat_y(jj, ii) = CompDom.Ucat_y(jj, mirror_loc);
            HaloDom.Pressure(jj, ii) = CompDom.Pressure(jj, mirror_loc);

            if ENABLE_DEBUGGING
                fprintf('DEBUG: \tGhost node (%d, %d) ', ii, jj);
                fprintf('mirrors inner node (%d, %d)\n', mirror_loc, jj);
            end
        end
    end

    %% Copy all values of ghost cells to the computational domain
    CompDom.Ucat_x = HaloDom.Ucat_x;
    CompDom.Ucat_y = HaloDom.Ucat_y;
    CompDom.Pressure = HaloDom.Pressure;

    %% Re-merge physical values to the inner domain
    for jj = ( 1+Nghost ):( N+Nghost )
        for ii = ( 1+Nghost ):( M+Nghost )
            CompDom.Ucat_x(jj, ii) = PhysDom.Ucat_x(jj-Nghost, ii-Nghost);
            CompDom.Ucat_y(jj, ii) = PhysDom.Ucat_y(jj-Nghost, ii-Nghost);
            CompDom.Pressure(jj, ii) = PhysDom.Pressure(jj-Nghost, ii-Nghost);
        end
    end

    fprintf('INFO: \tPeriodic boundary conditions are enabled!\n');
end

%% Return 'nan' to the corners of the computational domain