function [Ucont_x Ucont_y Ucat_x Ucat_y Ubcs_x Ubcs_y Pressure dx dy] = Init(M2,N2)
% This function initialize the driven cavity problem
L = 1;%pi
dx = L / (M2-3);
dy = L / (N2-3);

% Zero out everything
for i=1:M2
    for j = 1:N2
        
        x = (i-1.5) * dx;
        y = (j-1.5) * dy;
        
        Ucat_x(i,j) = 0;
        Ucat_y(i,j) = 0;

        saveCoordinates('coor_ucat.csv', x, y);
    end
end


% Setup contravariant component
for i = 1:M2+1
    for j = 1:N2
        
        x = (i-2) * dx;
        y = (j-1.5) * dy;                
        
        Ucont_x(i,j) =  0;

        saveCoordinates('coor_ucont_x.csv', x, y);
        
    end
end


for i = 1:M2
    for j = 1:N2+1
        
        x = (i-1.5) * dx;
        y = (j-2) * dy;        
        

        Ucont_y(i,j) = 0;

        saveCoordinates('coor_ucont_y.csv', x, y);
    end
end

% Set boundary condition
Ubcs_x(1:M2,1:N2)= 0;
Ubcs_y(1:M2,1:N2)= 0;

% Set the pressure
for i=1:M2
    for j = 1:N2
        
        x = (i-2) * dx;
        y = (j-2) * dy;       
    
        Pressure(i,j) = 0;

        saveCoordinates('coor_press.csv', x, y);
    end
end



