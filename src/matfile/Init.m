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
    dx = L / (M-1);
    dy = L / (N-1);

    %% Ghost variables' dimension in computational non-staggered grid
    M2 = M + (2*Nghost);
    N2 = N + (2*Nghost);

    %% Variables' dimension in staggered grid
    M3 = M + (2*Nghost) + 1;
    N3 = N + (2*Nghost) + 1;

    %% Allocate memory for the components
    Ucont_x = nan(N2, M2);
    Ucont_y = nan(N2, M2);

    Ucat_phys_x = nan(N, M);
    Ucat_phys_y = nan(N, M);
    Pressure_phys = nan(N, M);

    Ucat_halo_x = nan(N2, M2);
    Ucat_halo_y = nan(N2, M2);
    Pressure_halo = nan(N2, M2);

    Ubcs_x = nan(N2, M2);
    Ubcs_y = nan(N2, M2);

    %% Init face-centered velocity components
    %========= x-component =========
    for j = 1:N2
        for i = 1:M3
            
            x = (i -1 -Nghost -0.5) * dx;
            y = (j -1 -Nghost) * dy;                
            
            Ucont_x(j, i) =  U * sin( 2*pi*x ) * cos( 2*pi*y );

            if VISUAL_GRID
                %TAM_savecoorddinates('coord_face_centered_x.csv', x, y); % Old way, slow way
                coord_face_centered_x.x = [x coord_face_centered_x.x];
                coord_face_centered_x.y = [y coord_face_centered_x.y];
            end
        end
    end

    %========= y-component =========
    for j = 1:N3
        for i = 1:M2
            
            x = (i -1 -Nghost) * dx;
            y = (j -1 -Nghost -0.5) * dy;        
            
            Ucont_y(j, i) = -U * cos( 2*pi*x ) * sin( 2*pi*y );

            if VISUAL_GRID
                %TAM_savecoorddinates('coord_face_centered_y.csv', x, y); % Old way, slow way
                coord_face_centered_y.x = [x coord_face_centered_y.x];
                coord_face_centered_y.y = [y coord_face_centered_y.y];
            end
        end
    end

    %% Cell-centered velocity components
    [Ucat_comp_x, Ucat_comp_y] = Contra_To_Cart(Ucont_x, Ucont_y, M2, N2);

    Pressure_comp = nan(N2, M2);
    %% Init cell-centered components
    for j = 1:N2
        for i = 1:M2
            
            x = (i -1 -Nghost) * dx;
            y = (j -1 -Nghost) * dy;
            
            Pressure_comp(j, i) = -0.25 * U^2 * ( cos( 4*pi*x ) + cos( 4*pi*y ) );

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
                    x = (i -1 - Nghost) * dx;
                    y = (j -1 - Nghost) * dy;
                end

                %TAM_savecoorddinates('coord_ghost.csv', x, y); % Old way, slow way
                coord_ghost.x = [x coord_ghost.x];
                coord_ghost.y = [y coord_ghost.y];
            end
        end

        save('coordinates.mat', 'coord_cell_centered', 'coord_ghost', 'coord_face_centered_x', 'coord_face_centered_y');
    end

    %% Selecting the physical domain
    for jj = 1:N
        for ii = 1:M
            Ucat_phys_x(jj, ii) = Ucat_comp_x(jj+Nghost, ii+Nghost);
            Ucat_phys_y(jj, ii) = Ucat_comp_y(jj+Nghost, ii+Nghost);
            Pressure_phys(jj, ii) = Pressure_comp(jj+Nghost, ii+Nghost);
        end
    end

    %% Nesting the variables into structures
    PhysDom.Ucat_x = Ucat_phys_x;
    PhysDom.Ucat_y = Ucat_phys_y;
    PhysDom.Pressure = Pressure_phys;

    CompDom.Ucat_x = Ucat_comp_x;
    CompDom.Ucat_y = Ucat_comp_y;
    CompDom.Pressure = Pressure_comp;
    CompDom.Ucont_x = Ucont_x;
    CompDom.Ucont_y = Ucont_y;
    CompDom.Ubcs_x = Ubcs_x;
    CompDom.Ubcs_y = Ubcs_y;

    HaloDom.Ucat_x = Ucat_halo_x;
    HaloDom.Ucat_y = Ucat_halo_y;
    HaloDom.Pressure = Pressure_halo;

    %% Init fluxes
    FluxSum.Convective.Flux_x   = nan(N2, M2);
    FluxSum.Convective.Flux_y   = nan(N2, M2);
    FluxSum.Viscous.Flux_x      = nan(N2, M2);
    FluxSum.Viscous.Flux_y      = nan(N2, M2);
    FluxSum.P_Gradient.Flux_x   = nan(N2, M2);
    FluxSum.P_Gradient.Flux_y   = nan(N2, M2);
        
end