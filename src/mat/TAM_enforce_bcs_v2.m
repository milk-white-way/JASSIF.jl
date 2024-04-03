function [CompDomOut, HaloDom, Ucat_comp_x, Ucat_comp_y, Pressure_comp, Ucont_x, Ucont_y] = ...
    TAM_enforce_bcs_v2(M, N, M2, N2, Nghost, ...
                        PhysDom, CompDomIn, ...
                        iphys, iphye, jphys, jphye, ...
                        IS_PERIODIC, DEBUG)

    %% Version 2:
    %       + Is a function now
    %       + Less complicated periodic bcs procedure
    %       + Also includes the interpolation of the face-centered velocities

    %% Allocate local variables
    HaloDom.Ubcs_x = nan(M2, N2);
    HaloDom.Ubcs_y = nan(M2, N2);
    HaloDom.Pbcs   = nan(M2, N2);

    Ucont_x = CompDomIn.Ucont_x;
    Ucont_y = CompDomIn.Ucont_y;

    Ucat_comp_x   = nan(M2, N2);
    Ucat_comp_y   = nan(M2, N2);
    Pressure_comp = nan(M2, N2);

    Ucat_comp_x(iphys:iphye, jphys:jphye) = PhysDom.Ucat_x(1:M, 1:N);
    Ucat_comp_y(iphys:iphye, jphys:jphye) = PhysDom.Ucat_y(1:M, 1:N);
    Pressure_comp(iphys:iphye, jphys:jphye) = PhysDom.Pressure(1:M, 1:N);

    if IS_PERIODIC
        HaloDom.Ubcs_x(iphys:iphye,    1:Nghost) = Ucat_comp_x(iphys:iphye, jphye-Nghost+1:jphye); % South
        HaloDom.Ubcs_x(iphys:iphye,  jphye+1:N2) = Ucat_comp_x(iphys:iphye, jphys:jphys+Nghost-1); % North
        HaloDom.Ubcs_x(1:Nghost   , jphys:jphye) = Ucat_comp_x(iphye-Nghost+1:iphye, jphys:jphye); % West
        HaloDom.Ubcs_x(iphye+1:M2 , jphys:jphye) = Ucat_comp_x(iphys:iphys+Nghost-1, jphys:jphye); % East

        HaloDom.Ubcs_y(iphys:iphye,    1:Nghost) = Ucat_comp_y(iphys:iphye, jphye-Nghost+1:jphye); % South
        HaloDom.Ubcs_y(iphys:iphye,  jphye+1:N2) = Ucat_comp_y(iphys:iphye, jphys:jphys+Nghost-1); % North
        HaloDom.Ubcs_y(1:Nghost   , jphys:jphye) = Ucat_comp_y(iphye-Nghost+1:iphye, jphys:jphye); % West
        HaloDom.Ubcs_y(iphye+1:M2 , jphys:jphye) = Ucat_comp_y(iphys:iphys+Nghost-1, jphys:jphye); % East

        HaloDom.Pbcs(iphys:iphye,    1:Nghost) = Pressure_comp(iphys:iphye, jphye-Nghost+1:jphye); % South
        HaloDom.Pbcs(iphys:iphye,  jphye+1:N2) = Pressure_comp(iphys:iphye, jphys:jphys+Nghost-1); % North
        HaloDom.Pbcs(1:Nghost   , jphys:jphye) = Pressure_comp(iphye-Nghost+1:iphye, jphys:jphye); % West
        HaloDom.Pbcs(iphye+1:M2 , jphys:jphye) = Pressure_comp(iphys:iphys+Nghost-1, jphys:jphye); % East

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
            if ( ii <= Nghost) || (ii > (M + Nghost)) || (jj <= Nghost) || (jj > (N + Nghost) )
                Ucat_comp_x(ii, jj)   = HaloDom.Ubcs_x(ii, jj);
                Ucat_comp_y(ii, jj)   = HaloDom.Ubcs_y(ii, jj);
                Pressure_comp(ii, jj) = HaloDom.Pbcs(ii, jj);
            end
        end
    end

    % For Staggered grid
    Ucont_x(1, :)   = 0.5 * ( Ucat_comp_x(iphys, jphys:jphye) + Ucat_comp_x(iphys-1, jphys:jphye) ); % West
    Ucont_x(end, :) = 0.5 * ( Ucat_comp_x(iphye, jphys:jphye) + Ucat_comp_x(iphye+1, jphys:jphye) ); % East
    Ucont_y(:, 1)   = 0.5 * ( Ucat_comp_y(iphys:iphye, jphys) + Ucat_comp_y(iphys:iphye, jphys-1) ); % South
    Ucont_y(:, end) = 0.5 * ( Ucat_comp_y(iphys:iphye, jphye) + Ucat_comp_y(iphys:iphye, jphye+1) ); % North 

    % Ensemble the output structures
    CompDomOut.Ucat_x = Ucat_comp_x;
    CompDomOut.Ucat_y = Ucat_comp_y;
    CompDomOut.Pressure = Pressure_comp;
    CompDomOut.Ucont_x = Ucont_x;
    CompDomOut.Ucont_y = Ucont_y;

end