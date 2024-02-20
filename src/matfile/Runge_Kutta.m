function [FluxSumOut, U_im_x, U_im_y] = ...
    Runge_Kutta(CompDom, FluxSumIn, dU_old_x, dU_old_y, ...
                    M, N, M2, N2, Nghost, Re, dx, dy, dt, time, DEBUG)

    Ucont_x  = CompDom.Ucont_x; 
    Ucont_y  = CompDom.Ucont_y;
    Ucat_x   = CompDom.Ucat_x;
    Ucat_y   = CompDom.Ucat_y;
    Pressure = CompDom.Pressure;

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
            [Convective_Flux_x, Convective_Flux_y] = Convective_Flux(FluxSumIn, Ucont_x, Ucont_y, Ucat_x, Ucat_y, M2, N2, dx, dy);

            FluxSumOut.Convective.Flux_x = Convective_Flux_x;
            FluxSumOut.Convective.Flux_y = Convective_Flux_y;

            %% Calculate the viscous terms
            [Viscous_x, Viscous_y] = Viscous_Flux(FluxSumIn, Ucat_x, Ucat_y, M2, N2, dx, dy, Re, DEBUG);
            
            FluxSumOut.Viscous.Flux_x = Viscous_x;
            FluxSumOut.Viscous.Flux_y = Viscous_y;

            %% Calculate the pressure gradient terms
            [P_Gradient_x, P_Gradient_y] = Pressure_Gradient(FluxSumIn, Pressure, M2, N2, dx, dy);

            FluxSumOut.P_Gradient.Flux_x = P_Gradient_x;
            FluxSumOut.P_Gradient.Flux_y = P_Gradient_y;

            %% Calculate the total fluxes
            TotalFlux_x  = - FluxSumOut.Convective.Flux_x + FluxSumOut.Viscous.Flux_x - FluxSumOut.P_Gradient.Flux_x;
            TotalFlux_y  = - FluxSumOut.Convective.Flux_y + FluxSumOut.Viscous.Flux_y - FluxSumOut.P_Gradient.Flux_y;
        
            RHS_x = nan(N2, M2);
            RHS_y = nan(N2, M2);

            %% Then, fluxes are interpolated back to the face-centered to advance Ucont
            for jj = 2:N2-1
                for ii = 2:M2-1
                    RHS_x(jj, ii) = 0.5*(TotalFlux_x(jj, ii) + TotalFlux_x(jj, ii+1)) 
                    RHS_y(jj, ii) = 0.5*(TotalFlux_y(jj, ii) + TotalFlux_y(jj+1, ii))
                end
            end

            RHS_x = RHS_x - (1.5/dt) * (U_im_x  - Ucont_x) + (0.5/dt) * dU_old_x
            RHS_y = RHS_y - (1.5/dt) * (U_im_y  - Ucont_y) + (0.5/dt) * dU_old_y
        
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

end