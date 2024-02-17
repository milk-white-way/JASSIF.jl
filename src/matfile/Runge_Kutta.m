function [U_im_x, U_im_y] = ...
    Runge_Kutta(dU_old_x, dU_old_y, ...
                    CompDom, HaloDom, ...
                    M, N, M2, N2, M3, N3, ...
                    Re, dx, dy, dt, time)

Ucont_x  = CompDom.Ucont_x; 
Ucont_y  = CompDom.Ucont_y;
Ucat_x   = CompDom.Ucat_x;
Ucat_y   = CompDom.Ucat_y;
Pressure = CompDom.Pressure;

Ubcs_x   = HaloDom.Ubcs_x;
Ubcs_y   = HaloDom.Ubcs_y;

pseudo_t = 0;
alpha = [1/4; 1/3; 1/2; 1];

% Assign first guess
U_p_x = Ucont_x;
U_p_y = Ucont_y;

tol = 1e-8;
err = 1;

while (pseudo_t < 16 && err > tol)
     
    U_im_x = U_p_x;
    U_im_y = U_p_y;  
     
    for stage=1:4
        
        %% Fluxes computation first occurs at the cell-centered variables on non-staggered grid
        %[Convective_Flux_x, Convective_Flux_y] = Convective_Flux(Ucont_x, Ucont_y, Ucat_x, Ucat_y, dx, dy);
        %[Viscous_Flux_x, Viscous_Flux_y   ] = Viscous_Flux(Ucat_x, Ucat_y, dx, dy, Re);
        %[P_Gradient_x, P_Gradient_y] = Pressure_Gradient(Pressure, dx, dy);

        %% Fake fluxes
        FluxSum.Convective.Flux_x   = zeros(N2, M2);
        FluxSum.Convective.Flux_y   = zeros(N2, M2);
        FluxSum.Viscous.Flux_x      = zeros(N2, M2);
        FluxSum.Viscous.Flux_y      = zeros(N2, M2);
        FluxSum.P_Gradient.Flux_x   = zeros(N2, M2);
        FluxSum.P_Gradient.Flux_y   = zeros(N2, M2);

        TotalFlux_x  = - FluxSum.Convective.Flux_x + FluxSum.Viscous.Flux_x - FluxSum.P_Gradient.Flux_x; 
        TotalFlux_y  = - FluxSum.Convective.Flux_y + FluxSum.Viscous.Flux_y - FluxSum.P_Gradient.Flux_y;  
    
        %% Then, fluxes are interpolated back to the face-centered to advance Ucont
        for jj = 1:N2
            for ii = 1:M3
                RHS_x(jj, ii) = (TotalFlux_x(i,j) + TotalFlux_x(i+1,j)) / 2;  
            end
        end

        for i = 1:M
            for j= 1:N-1
                RHS_y(i,j) = (TotalFlux_y(i,j) + TotalFlux_y(i,j+1)) / 2;  
            end
        end

        RHS_x = zeros(N2, M3);
        RHS_y = zeros(N3, M2);

        RHS_x = RHS_x - (1.5/dt) * (U_im_x  - Ucont_x) + (0.5/dt) * dU_old_x;
        RHS_y = RHS_y - (1.5/dt) * (U_im_y  - Ucont_y) + (0.5/dt) * dU_old_y;
     
        U_im_x = U_p_x + alpha(stage) * dt * 0.4 * RHS_x;
        U_im_y = U_p_y + alpha(stage) * dt * 0.4 * RHS_y;         
                     
    end
         
    err = norm(U_p_x - U_im_x, inf); 
    U_p_x = U_im_x;
    U_p_y = U_im_y;    
     
    pseudo_t = pseudo_t+1;     
    
    fprintf('subitr = %d, Momentum convergence error = %8.6f \n', pseudo_t, err); 
end 
 
if (err > 1e-1)
    fprintf('time = %d, Momentum solver does not converge e = %8.6f \n', time, err);
end
 
 