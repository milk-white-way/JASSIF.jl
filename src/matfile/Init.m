function [Ucont_x, Ucont_y, Ucat_phy_x, Ucat_phy_y, Pressure_phy, Ucat_cal_x,  Ucat_cal_y, Pressure_cal, Ubcs_x, Ubcs_y, Pbcs, dx, dy] = Init(M, N, M2, N2, M3, N3)

L = 1; % Square domain
dx = L / M;
dy = L / N;

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
        
        %x = (i-1) * dx;
        %y = (j-1) * dy;
        
        Ucat_phy_x(i,j) = 0;
        Ucat_phy_y(i,j) = 0;
        Pressure_phy(i,j) = 0;

        %saveCoordinates('coor_cell_centered.csv', x, y);
    end
end

for i=1:M2
    for j = 1:N2

        %x = (i-1.5) * dx;
        %y = (j-1.5) * dy;

        if i == 1 || i == M2 || j == 1 || j == N2
            Ubcs_x(i,j) = 0;
            Ubcs_y(i,j) = 0;
            Pbcs(i,j) = 0;

            %saveCoordinates('coor_ghost.csv', x, y);
        end

    end
end

%% Init face-centered components
% x-component
for i = 1:M3
    for j = 1:N2
        
        %x = (i-2) * dx;
        %y = (j-1.5) * dy;                
        
        Ucont_x(i,j) =  0;

        %saveCoordinates('coor_ucont_x.csv', x, y);
        
    end
end

% y-component
for i = 1:M2
    for j = 1:N3
        
        %x = (i-1.5) * dx;
        %y = (j-2) * dy;        
        

        Ucont_y(i,j) = 0;

        %saveCoordinates('coor_ucont_y.csv', x, y);
    end
end

%% Create calculated domain
for i=1:M2
    for j = 1:N2

        %x = (i-1.5) * dx;
        %y = (j-1.5) * dy;

        if i == 1 || i == M2 || j == 1 || j == N2
            Ucat_cal_x(i,j) = Ubcs_x(i,j);
            Ucat_cal_y(i,j) = Ubcs_y(i,j);
            Pressure_cal(i,j) = Pbcs(i,j);
            %saveCoordinates('coor_ghost.csv', x, y);
        else 
            Ucat_cal_x(i,j) = Ucat_phy_x(i-1,j-1);
            Ucat_cal_y(i,j) = Ucat_phy_y(i-1,j-1);
            Pressure_cal(i,j) = Pressure_phy(i-1,j-1);
        end
    end
end