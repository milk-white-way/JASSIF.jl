function [Ucont_x, Ucont_y, Ucat_phy_x, Ucat_phy_y, Pressure_phy, Ucat_cal_x,  Ucat_cal_y, Pressure_cal, Ubcs_x, Ubcs_y, Pbcs, dx, dy] = Init(M, N, M2, N2, M3, N3, L, U, VISUAL_GRID)

    coord_cell_centered.x = [];
    coord_cell_centered.y = [];

    coord_ghost.x = [];
    coord_ghost.y = [];

    coord_face_centered_x.x = [];
    coord_face_centered_x.y = [];

    coord_face_centered_y.x = [];
    coord_face_centered_y.y = [];

    %% Initialization starts here
    dx = L / (M-1);
    dy = L / (N-1);

    %% Allocate memory for the components
    Ucont_x = nan(M3, N2);
    Ucont_y = nan(M2, N3);

    Ucat_phy_x = nan(M, N);
    Ucat_phy_y = nan(M, N);
    Pressure_phy = nan(M, N);

    Ubcs_x = nan(M2, N2);
    Ubcs_y = nan(M2, N2);
    Pbcs = nan(M2, N2);

    Ucat_cal_x = nan(M2, N2);
    Ucat_cal_y = nan(M2, N2);
    Pressure_cal = nan(M2, N2);

    %% Init cell-centered components, ghost cells included
    for i = 1:M
        for j = 1:N
            
            x = (i-1) * dx;
            y = (j-1) * dy;
            
            Ucat_phy_x(i,j) = U * sin( 2*pi*x ) * cos( 2*pi*y );
            Ucat_phy_y(i,j) = -U * cos( 2*pi*x ) * sin( 2*pi*y );
            Pressure_phy(i,j) = -0.25 * U^2 * ( cos( 4*pi*x ) + cos( 4*pi*y ) );

            if VISUAL_GRID
                TAM_savecoorddinates('coord_cell_centered.csv', x, y);
            end
            coord_cell_centered.x = [x coord_cell_centered.x];
            coord_cell_centered.y = [y coord_cell_centered.y];
        end
    end

    for i=1:M2
        for j = 1:N2

            x = (i-1-1) * dx;
            y = (j-1-1) * dy;

            if i == 1 || i == M2 || j == 1 || j == N2
                Ubcs_x(i,j) = U * sin( 2*pi*x ) * cos( 2*pi*y );
                Ubcs_y(i,j) = -U * cos( 2*pi*x ) * sin( 2*pi*y );
                Pbcs(i,j) = -0.25 * U^2 * ( cos( 4*pi*x ) + cos( 4*pi*y ) );

                if VISUAL_GRID
                    TAM_savecoorddinates('coord_ghost.csv', x, y);
                end
                coord_ghost.x = [x coord_ghost.x];
                coord_ghost.y = [y coord_ghost.y];
            end

        end
    end

    %% Init face-centered components
    % x-component
    for i = 1:M3
        for j = 1:N2
            
            x = (i-2-0.5) * dx;
            y = (j-2) * dy;                
            
            Ucont_x(i,j) =  U * sin( 2*pi*x ) * cos( 2*pi*y );

            if VISUAL_GRID
                TAM_savecoorddinates('coord_face_centered_x.csv', x, y);
            end
            coord_face_centered_x.x = [x coord_face_centered_x.x];
            coord_face_centered_x.y = [y coord_face_centered_x.y];
            
        end
    end

    % y-component
    for i = 1:M2
        for j = 1:N3
            
            x = (i-2) * dx;
            y = (j-2-0.5) * dy;        
            
            Ucont_y(i,j) = -U * cos( 2*pi*x ) * sin( 2*pi*y );

            if VISUAL_GRID
                TAM_savecoorddinates('coord_face_centered_y.csv', x, y);
            end
            coord_face_centered_y.x = [x coord_face_centered_y.x];
            coord_face_centered_y.y = [y coord_face_centered_y.y];
        end
    end

    %% Create calculated domain by merging physical and ghost cells
    for i=1:M2
        for j = 1:N2

            if i == 1 || i == M2 || j == 1 || j == N2
                Ucat_cal_x(i,j) = Ubcs_x(i,j);
                Ucat_cal_y(i,j) = Ubcs_y(i,j);
                Pressure_cal(i,j) = Pbcs(i,j);
            else 
                Ucat_cal_x(i,j) = Ucat_phy_x(i-1,j-1);
                Ucat_cal_y(i,j) = Ucat_phy_y(i-1,j-1);
                Pressure_cal(i,j) = Pressure_phy(i-1,j-1);
            end
        end
    end

    save('coordinates.mat', 'coord_cell_centered', 'coord_ghost', 'coord_face_centered_x', 'coord_face_centered_y');

end