function [FluxSumOut, U_im_x, U_im_y] = ...
    Runge_Kutta(CompDom, FluxSumIn, dU_old_x, dU_old_y, ...
                    M, N, M2, N2, M3, N3, Nghost, Re, dx, dy, dt, time, ... 
                    iphys, iphye, jphys, jphye, DEBUG)

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
            [Convective_Flux_x, Convective_Flux_y] = Convective_Flux(FluxSumIn, Ucont_x, Ucont_y, Ucat_x, Ucat_y, M, N, Nghost, dx, dy, iphys, iphye, jphys, jphye, DEBUG);

            FluxSumOut.Convective.Flux_x = Convective_Flux_x(iphys:iphye, jphys:jphye);
            FluxSumOut.Convective.Flux_y = Convective_Flux_y(iphys:iphye, jphys:jphye);

            %% Calculate the viscous terms
            [Viscous_x, Viscous_y] = Viscous_Flux(FluxSumIn, Ucat_x, Ucat_y, M, N, Nghost, dx, dy, Re, DEBUG);
            
            FluxSumOut.Viscous.Flux_x = Viscous_x(iphys:iphye, jphys:jphye);
            FluxSumOut.Viscous.Flux_y = Viscous_y(iphys:iphye, jphys:jphye);

            %% Calculate the pressure gradient terms
            [P_Gradient_x, P_Gradient_y] = Pressure_Gradient(FluxSumIn, Pressure, M, N, Nghost, dx, dy);

            FluxSumOut.P_Gradient.Flux_x = P_Gradient_x(iphys:iphye, jphys:jphye);
            FluxSumOut.P_Gradient.Flux_y = P_Gradient_y(iphys:iphye, jphys:jphye);

            %% Calculate the total fluxes
            TotalFlux_x  = - Convective_Flux_x + Viscous_x - P_Gradient_x;
            TotalFlux_y  = - Convective_Flux_y + Viscous_y - P_Gradient_y;
            
            %% DEBATE: Apply boundary conditions to the fluxes
            TotalFlux_x(1:Nghost   , iphys:iphye) = TotalFlux_x(jphye-Nghost+1:jphye, iphys:iphye); % South
            TotalFlux_x(jphye+1:N2 , iphys:iphye) = TotalFlux_x(jphys:jphys+Nghost-1, iphys:iphye); % North
            TotalFlux_x(jphys:jphye,    1:Nghost) = TotalFlux_x(jphys:jphye, iphye-Nghost+1:iphye); % West
            TotalFlux_x(jphys:jphye,  iphye+1:M2) = TotalFlux_x(jphys:jphye, iphys:iphys+Nghost-1); % East


            TotalFlux_y(1:Nghost   , iphys:iphye) = TotalFlux_y(jphye-Nghost+1:jphye, iphys:iphye); % South
            TotalFlux_y(jphye+1:N2 , iphys:iphye) = TotalFlux_y(jphys:jphys+Nghost-1, iphys:iphye); % North
            TotalFlux_y(jphys:jphye,    1:Nghost) = TotalFlux_y(jphys:jphye, iphye-Nghost+1:iphye); % West
            TotalFlux_y(jphys:jphye,  iphye+1:M2) = TotalFlux_y(jphys:jphye, iphys:iphys+Nghost-1); % East

            %% Then, fluxes are interpolated back to the face-centered to advance Ucont
            % 2-Step Process
            RHS_comp_x = nan(N2, M2-1);
            RHS_comp_y = nan(N2-1, M2);

            for ii = 1:M2-1
                for jj = 1:N2
                    RHS_comp_x(jj, ii) = 0.5*(TotalFlux_x(jj, ii) + TotalFlux_x(jj, ii+1));
                end
            end

            for ii = 1:M2
                for jj = 1:N2-1
                    RHS_comp_y(jj, ii) = 0.5*(TotalFlux_y(jj, ii) + TotalFlux_y(jj+1, ii));
                end
            end

            RHS_phys_x = RHS_comp_x(jphys:jphye, iphys-1:iphye);
            RHS_phys_y = RHS_comp_y(jphys-1:jphye, iphys:iphye);

            RHS_phys_x = RHS_phys_x - (1.5/dt) * (U_im_x  - Ucont_x) + (0.5/dt) * dU_old_x;
            RHS_phys_y = RHS_phys_y - (1.5/dt) * (U_im_y  - Ucont_y) + (0.5/dt) * dU_old_y;
        
            U_im_x = U_p_x + alpha(stage) * dt * 0.4 * RHS_phys_x;
            U_im_y = U_p_y + alpha(stage) * dt * 0.4 * RHS_phys_y;         
                        
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