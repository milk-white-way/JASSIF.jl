function [PhysDom, CompDom, HaloDom, FluxSum, M2, N2, dx, dy] = ...
    Init(M, N, Nghost, L, U, VISUAL_GRID)
    if VISUAL_GRID
        coord_cell_centered.x = [];
        coord_cell_centered.y = [];

        coord_ghost.x = [];
        coord_ghost.y = [];

        coord_face_centered_x.x = [];
        coord_face_centered_x.y = [];

        coord_face_centered_y.x = [];
        coord_face_centered_y.y = [];
    end

    %% Initialization starts here
    dx = L / M;
    dy = L / N;

    %% Ghost variables' dimension in computational non-staggered grid
    M2 = M + (2*Nghost);
    N2 = N + (2*Nghost);

    %% Variables' dimension in staggered grid
    %M3 = M + (2*Nghost) + 1;
    %N3 = N + (2*Nghost) + 1;
    M3 = M+1;
    N3 = N+1;

    %% Allocate memory for the components
    Ucont_x = nan(N, M3);
    Ucont_y = nan(N3, M);

    Pressure_phys = nan(N, M);

    Ucat_halo_x = nan(N2, M2);
    Ucat_halo_y = nan(N2, M2);
    Pressure_halo = nan(N2, M2);

    %% Init face-centered velocity components
    %========= x-component =========
    for i = 1:M3 
        for j = 1:N
            
            x = (i -1) * dx;
            y = (j -1 + 0.5) * dy;                
            
            Ucont_x(j, i) =  U * sin( 2*pi*x ) * cos( 2*pi*y );

            if VISUAL_GRID
                %TAM_savecoorddinates('coord_face_centered_x.csv', x, y); % Old way, slow way
                coord_face_centered_x.x = [x coord_face_centered_x.x];
                coord_face_centered_x.y = [y coord_face_centered_x.y];
            end
        end
    end

    %========= y-component =========
    for i = 1:M
        for j = 1:N3
            
            x = (i -1 +0.5) * dx;
            y = (j -1) * dy;        
            
            Ucont_y(j, i) = -U * cos( 2*pi*x ) * sin( 2*pi*y );

            if VISUAL_GRID
                %TAM_savecoorddinates('coord_face_centered_y.csv', x, y); % Old way, slow way
                coord_face_centered_y.x = [x coord_face_centered_y.x];
                coord_face_centered_y.y = [y coord_face_centered_y.y];
            end
        end
    end

    %% Convert to get the physical domain
    [Ucat_phys_x, Ucat_phys_y] = Contra_To_Cart(Ucont_x, Ucont_y, M, N);

    %% Init cell-centered components
    for i = 1:M
        for j = 1:N
            
            x = (i -1 +0.5) * dx;
            y = (j -1 +0.5) * dy;
            
            Pressure_phys(j, i) = -0.25 * U^2 * ( cos( 4*pi*x ) + cos( 4*pi*y ) );

            if VISUAL_GRID
                %TAM_savecoorddinates('coord_cell_centered.csv', x, y); % Old way, slow way
                coord_cell_centered.x = [x coord_cell_centered.x];
                coord_cell_centered.y = [y coord_cell_centered.y];
            end
        end
    end

    if VISUAL_GRID
        %% Init the ghost cells here only neccessitate debugging
        for j = 1:N2
            for i = 1:M2
                if i <= Nghost || i > (M2-Nghost) || j <= Nghost || j > (N2-Nghost)
                    x = (i -1 - Nghost +0.5) * dx;
                    y = (j -1 - Nghost +0.5) * dy;
                end

                %TAM_savecoorddinates('coord_ghost.csv', x, y); % Old way, slow way
                coord_ghost.x = [x coord_ghost.x];
                coord_ghost.y = [y coord_ghost.y];
            end
        end

        save('coordinates.mat', 'coord_cell_centered', 'coord_ghost', 'coord_face_centered_x', 'coord_face_centered_y');
    end

    %% Nesting the variables into structures
    PhysDom.Ucat_x = Ucat_phys_x;
    PhysDom.Ucat_y = Ucat_phys_y;
    PhysDom.Pressure = Pressure_phys;

    HaloDom.Ucat_x = Ucat_halo_x;
    HaloDom.Ucat_y = Ucat_halo_y;
    HaloDom.Pressure = Pressure_halo;

    [Ucat_comp_x, Ucat_comp_y, Pressure_comp] = TAM_ensemble_comp(PhysDom, HaloDom, M, N, Nghost);

    CompDom.Ucat_x = Ucat_comp_x;
    CompDom.Ucat_y = Ucat_comp_y;
    CompDom.Pressure = Pressure_comp;
    CompDom.Ucont_x = Ucont_x;
    CompDom.Ucont_y = Ucont_y;

    %% Init fluxes
    FluxSum.Convective.Flux_x   = nan(N2, M2);
    FluxSum.Convective.Flux_y   = nan(N2, M2);
    FluxSum.Viscous.Flux_x      = nan(N2, M2);
    FluxSum.Viscous.Flux_y      = nan(N2, M2);
    FluxSum.P_Gradient.Flux_x   = nan(N2, M2);
    FluxSum.P_Gradient.Flux_y   = nan(N2, M2);
        
end