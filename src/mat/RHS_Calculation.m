function [RHS_x, RHS_y] = ...
    RHS_Calculation(Ucont_x, Ucont_y, ...
                    Ucat_x, Ucat_y, ...
                    M, N, M2, N2, M3, N3, ...
                    Pressure, Re,dx, dy)

    %[Convective_Flux_x, Convective_Flux_y] = Convective_Flux(Ucont_x, Ucont_y, Ucat_x, Ucat_y, dx, dy);
    %[Viscous_Flux_x, Viscous_Flux_y   ] = Viscous_Flux(Ucat_x, Ucat_y, dx, dy, Re);
    %[P_Gradient_x, P_Gradient_y] = Pressure_Gradient(Pressure, dx, dy);

    %% Fake fluxes
    FluxFld.Convective_Flux_x = zeros(N2, N2);
    FluxFld.Convective_Flux_y = zeros(N2, N2);
    FluxFld.Viscous_Flux_x = zeros(N2, N2);
    FluxFld.Viscous_Flux_y = zeros(N2, N2);
    FluxFld.P_Gradient_x = zeros(N2, N2);
    FluxFld.P_Gradient_y = zeros(N2, N2);

    F_x  = - FluxFld.Convective_Flux_x + FluxFld.Viscous_Flux_x - FluxFld.P_Gradient_x;       
    F_y  = - FluxFld.Convective_Flux_y + FluxFld.Viscous_Flux_y - FluxFld.P_Gradient_y;       
    
    % Interpolation back to the half node formulation to advance Ucont
    for i = 1:M-1
        for j = 1:N
            RHS_x(i,j) = (F_x(i,j) + F_x(i+1,j)) / 2;  
        end
    end
    
    % Zero out the boundaries
    for j=1:N
        RHS_x(M,j) = 0;
        RHS_x(M-1,j) = 0;
        RHS_x(1,j) = 0;
    end

    for i = 1:M
        RHS_x(i,1) = 0;
        RHS_x(i,N) = 0;
    end

    for i = 1:M
        for j= 1:N-1
            RHS_y(i,j) = (F_y(i,j) + F_y(i,j+1)) / 2;  
        end
    end
    
    % Zero out the boundary
    for i=1:M
        RHS_y(i,N-1) = 0;
        RHS_y(i,1) = 0;
        RHS_y(i,N) = 0;
    end
    
    for j=1:N
        RHS_y(1,j) = 0;
        RHS_y(M,j) = 0;
    end
end 