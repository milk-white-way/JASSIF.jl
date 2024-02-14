function [Ucont_x Ucont_y Ucat_x Ucat_y Ubcs_x Ubcs_y Pressure dx dy] = Init_Taylor_Green(M2,N2)
% This function initialize the driven cavity problem
L = pi;
dx = pi / (M2-3);
dy = pi / (N2-3);

% Zero out everything
for i=1:M2
    for j = 1:N2
        
        x = (i-2) * dx;
        y = (j-2) * dy;
        
    Ucat_x(i,j) = -cos(x) * sin(y);
    Ucat_y(i,j) = sin(x) * cos(y);
    end
end


% Setup contravariant component
for i = 1:M2
    for j = 1:N2
        
        x = (i-2+0.5) * dx;
        y = (j-2) * dy;                
        
        Ucont_x(i,j) =  -cos(x) * sin(y);
        
    end
end


for i = 1:M2
    for j = 1:N2
        
        x = (i-2) * dx;
        y = (j-2+0.5) * dy;        
        

        Ucont_y(i,j) = sin(x) * cos(y);
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
    
    Pressure(i,j) = -(cos(2*x) + cos(2*y))/4;
    end
end



