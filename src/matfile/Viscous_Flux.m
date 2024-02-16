function [Viscous_x, Viscous_y] = Viscous_Flux(Ucat_x, Ucat_y, dx, dy, Re)
% Re is Reynolds number
M = length(Ucat_x(:,1));
N = length(Ucat_x(1,:));

% Take the central differencing of this thing
for i = 2 :M-1
    for j = 2:N-1
        
        % ------------- x direction -------------------
        % E-W
        ue = Ucat_x(i-1,j);
        up = Ucat_x(i,j);
        uw = Ucat_x(i+1,j);
        
        % N-S
        un = Ucat_x(i,j-1);
        up = Ucat_x(i,j);
        us = Ucat_x(i,j+1);
        
        Viscous_x(i,j) = (uw - 2*up + ue)./(dx^2) + (us - 2*up + un)./(dy^2);
        % ------------- y direction --------------------
        % E-W
        ve = Ucat_y(i-1,j);
        vp = Ucat_y(i,j);
        vw = Ucat_y(i+1,j);
        
        % N-S
        vn = Ucat_y(i,j-1);
        vp = Ucat_y(i,j);
        vs = Ucat_y(i,j+1);
        
        Viscous_y(i,j) = (vw - 2*vp + ve)./(dx^2) + (vs - 2*vp + vn)./(dy^2);
    end
end

% At the boundary viscous flux = 0
for j=1:N
    Viscous_x(1,j) = 0;
    Viscous_x(M,j) = 0;
    
    Viscous_y(1,j) = 0;    
    Viscous_y(M,j) = 0;
end

for i = 1:M
    Viscous_x(i,1) = 0;
    Viscous_x(i,N) = 0;
    
    Viscous_y(i,1) = 0;
    Viscous_y(i,N) = 0;
end

% Scale with Reynolds number
Viscous_x = Viscous_x / Re;
Viscous_y = Viscous_y / Re;