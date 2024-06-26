function [PhysDom, CompDom, M2, N2, M3, N3, iphys, iphye, jphys, jphye, dx, dy] = ...
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

    %% Quality of life here
    %==================================
    iphys = 1+Nghost; iphye = M+Nghost;
    jphys = 1+Nghost; jphye = N+Nghost;
    %==================================

    %% Allocate memory for the components
    Ucont_x = nan(M3, N);
    Ucont_y = nan(M, N3);

    Pressure_phys = nan(M, N);

    %% Init face-centered velocity components
    %========= x-component =========
    for i = 1:M3 
        for j = 1:N
            
            x = (i -1) * dx;
            y = (j -1 + 0.5) * dy;                
            
            Ucont_x(i, j) =  U * sin( 2*pi*x ) * cos( 2*pi*y );

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
            
            Ucont_y(i, j) = -U * cos( 2*pi*x ) * sin( 2*pi*y );

            if VISUAL_GRID
                %TAM_savecoorddinates('coord_face_centered_y.csv', x, y); % Old way, slow way
                coord_face_centered_y.x = [x coord_face_centered_y.x];
                coord_face_centered_y.y = [y coord_face_centered_y.y];
            end
        end
    end

    %% Convert to get the physical domain
    [Ucat_phys_x, Ucat_phys_y] = Contra_To_Cart(M, N, Ucont_x, Ucont_y);

    %% Init cell-centered components
    for i = 1:M
        for j = 1:N
            
            x = (i -1 +0.5) * dx;
            y = (j -1 +0.5) * dy;
            
            Pressure_phys(i, j) = 0.25 * U^2 * ( cos( 4*pi*x ) + cos( 4*pi*y ) );

            if VISUAL_GRID
                %TAM_savecoorddinates('coord_cell_centered.csv', x, y); % Old way, slow way
                coord_cell_centered.x = [x coord_cell_centered.x];
                coord_cell_centered.y = [y coord_cell_centered.y];
            end
        end
    end

    if VISUAL_GRID
        for i = 1:M2
            for j = 1:N2
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

    CompDom.Ucont_x = Ucont_x;
    CompDom.Ucont_y = Ucont_y;
 
end