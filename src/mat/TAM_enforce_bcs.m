fprintf('\nINFO: \tEnforcing boundary conditions... \n');

if ENABLE_BC_PERIODIC
    %% North and South boundaries:
    % North wall is the same as the South wall
    %HaloDom.Ubcs_x(1+Nghost:N+Nghost, Nghost+1) = PhysDom.Ucat_x(:, M);
    %HaloDom.Ubcs_y(1+Nghost:N+Nghost, Nghost+1) = PhysDom.Ucat_y(:, M);
    %HaloDom.Ubcs_x(1+Nghost:N+Nghost, M+Nghost) = PhysDom.Ucat_x(:, 1);
    %HaloDom.Ubcs_y(1+Nghost:N+Nghost, M+Nghost) = PhysDom.Ucat_y(:, 1);

    % Halo region
    for ii = iphys:iphye
        for jj = 1:Nghost 
            mirror_dis = Nghost - jj;
            mirror_loc = N + Nghost - mirror_dis;
            %mirror_loc = Nghost + 1 + mirror_dis;

            HaloDom.Ubcs_x(jj, ii) = CompDom.Ucat_x(mirror_loc, ii);
            HaloDom.Ubcs_y(jj, ii) = CompDom.Ucat_y(mirror_loc, ii);
            HaloDom.Pbcs(jj, ii) = CompDom.Pressure(mirror_loc, ii);

            if ENABLE_DEBUGGING
                fprintf('DEBUG: \tGhost node (%d, %d) ', ii, jj);
                fprintf('mirrored inner node (%d, %d)\n', ii, mirror_loc);
            end
        end

        for jj = (N+Nghost+1):N2
            mirror_dis = jj - (N + 1 + Nghost);
            mirror_loc = 1 + Nghost + mirror_dis;
            %mirror_loc = N + Nghost + mirror_dis;

            HaloDom.Ubcs_x(jj, ii) = CompDom.Ucat_x(mirror_loc, ii);
            HaloDom.Ubcs_y(jj, ii) = CompDom.Ucat_y(mirror_loc, ii);
            HaloDom.Pbcs(jj, ii) = CompDom.Pressure(mirror_loc, ii);

            if ENABLE_DEBUGGING
                fprintf('DEBUG: \tGhost node (%d, %d) ', ii, jj);
                fprintf('mirrored inner node (%d, %d)\n', ii, mirror_loc);
            end
        end
    end

    %% East and West boundaries:
    % East wall is the same as the West wall
    %HaloDom.Ubcs_x(Nghost+1, 1+Nghost:M+Nghost) = PhysDom.Ucat_x(N, :);
    %HaloDom.Ubcs_y(Nghost+1, 1+Nghost:M+Nghost) = PhysDom.Ucat_y(N, :);
    %HaloDom.Ubcs_x(N+Nghost, 1+Nghost:M+Nghost) = PhysDom.Ucat_x(1, :);
    %HaloDom.Ubcs_y(N+Nghost, 1+Nghost:M+Nghost) = PhysDom.Ucat_y(1, :);

    for jj = jphys:jphye
        for ii = 1:Nghost
            mirror_dis = Nghost - ii;
            mirror_loc = M + Nghost - mirror_dis;

            HaloDom.Ubcs_x(jj, ii) = CompDom.Ucat_x(jj, mirror_loc);
            HaloDom.Ubcs_y(jj, ii) = CompDom.Ucat_y(jj, mirror_loc);
            HaloDom.Pbcs(jj, ii) = CompDom.Pressure(jj, mirror_loc);

            if ENABLE_DEBUGGING
                fprintf('DEBUG: \tGhost node (%d, %d) ', ii, jj);
                fprintf('mirrored inner node (%d, %d)\n', mirror_loc, jj);
            end
        end

        for ii = (M+Nghost+1):M2
            mirror_dis = ii - (M + 1 + Nghost);
            mirror_loc = 1 + Nghost + mirror_dis;

            HaloDom.Ubcs_x(jj, ii) = CompDom.Ucat_x(jj, mirror_loc);
            HaloDom.Ubcs_y(jj, ii) = CompDom.Ucat_y(jj, mirror_loc);
            HaloDom.Pbcs(jj, ii) = CompDom.Pressure(jj, mirror_loc);

            if ENABLE_DEBUGGING
                fprintf('DEBUG: \tGhost node (%d, %d) ', ii, jj);
                fprintf('mirrored inner node (%d, %d)\n', mirror_loc, jj);
            end
        end
    end

    fprintf('INFO: \tPERIODIC BOUNDARY CONDITIONS ENFORCED\n');
end

[Ucat_comp_x, Ucat_comp_y, Pressure_comp] = TAM_ensemble_comp(PhysDom, HaloDom, M, N, Nghost);

CompDom.Ucat_x = Ucat_comp_x;
CompDom.Ucat_y = Ucat_comp_y;
CompDom.Pressure = Pressure_comp;