function [CompDomOut, HaloDom, Ucat_comp_x, Ucat_comp_y, Pressure_comp, Ucont_x, Ucont_y] = ...
    TAM_enforce_bcs_v2(PhysDom, CompDomIn, M, N, M2, N2, Nghost, iphys, iphye, jphys, jphye, IS_PERIODIC, DEBUG)

    %% Version 2:
    %       + Is a function now
    %       + Less complicated periodic bcs procedure
    %       + Also includes the interpolation of the face-centered velocities

    %% Allocate local variables
    HaloDom.Ubcs_x = nan(N2,M2);
    HaloDom.Ubcs_y = nan(N2,M2);
    HaloDom.Pbcs   = nan(N2,M2);

    Ucont_x = CompDomIn.Ucont_x;
    Ucont_y = CompDomIn.Ucont_y;

    Ucat_comp_x   = nan(N2,M2);
    Ucat_comp_y   = nan(N2,M2);
    Pressure_comp = nan(N2,M2);

    Ucat_comp_x(jphys:jphye, iphys:iphye) = PhysDom.Ucat_x(1:N, 1:M);
    Ucat_comp_y(jphys:jphye, iphys:iphye) = PhysDom.Ucat_y(1:N, 1:M);
    Pressure_comp(jphys:jphye, iphys:iphye) = PhysDom.Pressure(1:N, 1:M);

    if IS_PERIODIC
        HaloDom.Ubcs_x(1:Nghost   , iphys:iphye) = Ucat_comp_x(jphye-Nghost+1:jphye, iphys:iphye); % South
        HaloDom.Ubcs_x(jphye+1:N2 , iphys:iphye) = Ucat_comp_x(jphys:jphys+Nghost-1, iphys:iphye); % North
        HaloDom.Ubcs_x(jphys:jphye,    1:Nghost) = Ucat_comp_x(jphys:jphye, iphye-Nghost+1:iphye); % West
        HaloDom.Ubcs_x(jphys:jphye,  iphye+1:M2) = Ucat_comp_x(jphys:jphye, iphys:iphys+Nghost-1); % East

        HaloDom.Ubcs_y(1:Nghost   , iphys:iphye) = Ucat_comp_y(jphye-Nghost+1:jphye, iphys:iphye); % South
        HaloDom.Ubcs_y(jphye+1:N2 , iphys:iphye) = Ucat_comp_y(jphys:jphys+Nghost-1, iphys:iphye); % North
        HaloDom.Ubcs_y(jphys:jphye,    1:Nghost) = Ucat_comp_y(jphys:jphye, iphye-Nghost+1:iphye); % West
        HaloDom.Ubcs_y(jphys:jphye,  iphye+1:M2) = Ucat_comp_y(jphys:jphye, iphys:iphys+Nghost-1); % East

        HaloDom.Pbcs(1:Nghost   , iphys:iphye)   = Pressure_comp(jphye-Nghost+1:jphye, iphys:iphye); % South
        HaloDom.Pbcs(jphye+1:N2 , iphys:iphye)   = Pressure_comp(jphys:jphys+Nghost-1, iphys:iphye); % North
        HaloDom.Pbcs(jphys:jphye,    1:Nghost)   = Pressure_comp(jphys:jphye, iphye-Nghost+1:iphye); % West
        HaloDom.Pbcs(jphys:jphye,  iphye+1:M2)   = Pressure_comp(jphys:jphye, iphys:iphys+Nghost-1); % East

        if DEBUG
            fprintf('INFO: \tPERIODIC BOUNDARY CONDITIONS ENFORCED SUCCESSFULLY!\n');
            fprintf('INFO: \tComputational domain will be reconstructed.\n');
        end
    else
        warning('INFO: \tWALL BOUNDARY CONDITIONS NOT IMPLEMENTED YET!\n');
    end

    %% Reconstruct the computational domain with the ghost cells
    % For Non-staggered grid
    for ii = 1:M2
        for jj = 1:N2
            if (ii <= Nghost) || (ii > (M + Nghost)) || (jj <= Nghost) || (jj > (N + Nghost))
                Ucat_comp_x(jj, ii)   = HaloDom.Ubcs_x(jj, ii);
                Ucat_comp_y(jj, ii)   = HaloDom.Ubcs_y(jj, ii);
                Pressure_comp(jj, ii) = HaloDom.Pbcs(jj, ii);
            end
        end
    end

    % For Staggered grid
    Ucont_x(:, 1)   = 0.5 * ( Ucat_comp_x(jphys:jphye, iphys) + Ucat_comp_x(jphys:jphye, iphys-1) ); % West
    Ucont_x(:, end) = 0.5 * ( Ucat_comp_x(jphys:jphye, iphye) + Ucat_comp_x(jphys:jphye, iphye+1) ); % East
    Ucont_y(1, :)   = 0.5 * ( Ucat_comp_y(jphys, iphys:iphye) + Ucat_comp_y(jphys-1, iphys:iphye) ); % South
    Ucont_y(end, :) = 0.5 * ( Ucat_comp_y(jphye, iphys:iphye) + Ucat_comp_y(jphye+1, iphys:iphye) ); % North 

    % Ensemble the output structures
    CompDomOut.Ucat_x = Ucat_comp_x;
    CompDomOut.Ucat_y = Ucat_comp_y;
    CompDomOut.Pressure = Pressure_comp;
    CompDomOut.Ucont_x = Ucont_x;
    CompDomOut.Ucont_y = Ucont_y;

end