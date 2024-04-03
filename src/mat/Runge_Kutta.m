function [FluxSumOut, Ucont_im_x, Ucont_im_y] = ...
    Runge_Kutta(M, N, M2, N2, M3, N3, Nghost, ...
                CompDom, dU_old_x, dU_old_y, ...
                iphys, iphye, jphys, jphye, ...
                Re, dx, dy, dt, time, IS_PERIODIC, DEBUG)

    Ucont_prev_x = CompDom.Ucont_x;
    Ucont_prev_y = CompDom.Ucont_y;

    pseudo_t = 0;
    alpha = [1/4; 1/3; 1/2; 1];

    % Assign initial guesses
    Ucont_ig_x = Ucont_prev_x;
    Ucont_ig_y = Ucont_prev_y;

    Ucat_ig_x   = CompDom.Ucat_x;
    Ucat_ig_y   = CompDom.Ucat_y;
    Pressure_ig = CompDom.Pressure;

    tol = 1e-12;
    err = 1;

    while (pseudo_t < 100 && err > tol)
        
        Ucont_im_x = Ucont_ig_x;
        Ucont_im_y = Ucont_ig_y; 

        Ucat_im_x = Ucat_ig_x;
        Ucat_im_y = Ucat_ig_y;
        Pressure_im = Pressure_ig;
        
        for stage=1:4
            
            %% Fluxes computation first occurs at the cell-centered variables on non-staggered grid
            [Convective_x, Convective_y] = TAM_Convective_Flux(M, N, M2, N2, M3, N3, Nghost, Ucont_im_x, Ucont_im_y, Ucat_im_x, Ucat_im_y, iphys, iphye, jphys, jphye, dx, dy, DEBUG);

            FluxSumOut.Convective.Flux_x = Convective_x(iphys:iphye, jphys:jphye);
            FluxSumOut.Convective.Flux_y = Convective_y(iphys:iphye, jphys:jphye);

            % Calculate the viscous terms
            [Viscous_x, Viscous_y] = Viscous_Flux(M2, N2, Ucat_im_x, Ucat_im_y, iphys, iphye, jphys, jphye, dx, dy, Re, DEBUG);
            
            FluxSumOut.Viscous.Flux_x = Viscous_x(iphys:iphye, jphys:jphye);
            FluxSumOut.Viscous.Flux_y = Viscous_y(iphys:iphye, jphys:jphye);

            % Calculate the pressure gradient terms
            [P_Gradient_x, P_Gradient_y] = Pressure_Gradient(M2, N2, Pressure_im, iphys, iphye, jphys, jphye, dx, dy, DEBUG);

            FluxSumOut.P_Gradient.Flux_x = P_Gradient_x(iphys:iphye, jphys:jphye);
            FluxSumOut.P_Gradient.Flux_y = P_Gradient_y(iphys:iphye, jphys:jphye);

            % Calculate the total fluxes
            TotalFlux_x  = - FluxSumOut.Convective.Flux_x + FluxSumOut.Viscous.Flux_x - FluxSumOut.P_Gradient.Flux_x;
            TotalFlux_y  = - FluxSumOut.Convective.Flux_y + FluxSumOut.Viscous.Flux_y - FluxSumOut.P_Gradient.Flux_y;
            
            %% Then, fluxes are interpolated back to the face-centered to advance intermediate Ucont components from their initial guesses
            % Right hand side
            RHS_x = zeros(M3, N);
            RHS_y = zeros(M, N3);

            for ii = 2:M
                for jj = 1:N
                    RHS_x(ii, jj) = 0.5*( TotalFlux_x(ii-1, jj) + TotalFlux_x(ii, jj) );
                end
            end

            for ii = 1:M
                for jj = 2:N
                    RHS_y(ii, jj) = 0.5*( TotalFlux_y(ii, jj-1) + TotalFlux_y(ii, jj) );
                end
            end

            RHS_x = RHS_x - (1.5/dt) * (Ucont_im_x  - Ucont_prev_x) + (0.5/dt) * dU_old_x;
            RHS_y = RHS_y - (1.5/dt) * (Ucont_im_y  - Ucont_prev_y) + (0.5/dt) * dU_old_y;
        
            Ucont_im_x = Ucont_ig_x + alpha(stage) * dt * 0.4 * RHS_x;
            Ucont_im_y = Ucont_ig_y + alpha(stage) * dt * 0.4 * RHS_y;

            %% Finally, Update the boundary conditions at the next pseudo time step
            % Convert the intermediate Ucont components to intermediate Ucat components
            [Ucat_in_x, Ucat_in_y] = Contra_To_Cart(M, N, Ucont_im_x, Ucont_im_y);
            % Update the boundary conditions on the intermediate Ucat components
            PseudoPhysDom.Ucat_x = Ucat_in_x;
            PseudoPhysDom.Ucat_y = Ucat_in_y;
            PseudoPhysDom.Pressure = Pressure_im(iphys:iphye, jphys:jphye);
            PseudoCompDom.Ucont_x = Ucont_im_x;
            PseudoCompDom.Ucont_y = Ucont_im_y;

            [~, ~, Ucat_im_x, Ucat_im_y, Pressure_im, Ucont_im_x, Ucont_im_y] = TAM_enforce_bcs_v2(M, N, M2, N2, Nghost, PseudoPhysDom, PseudoCompDom, iphys, iphye, jphys, jphye, IS_PERIODIC, DEBUG);

        end

        erx = norm(Ucont_ig_x - Ucont_im_x, inf);
        ery = norm(Ucont_ig_y - Ucont_im_y, inf);
        err = max(erx, ery);

        Ucont_ig_x = Ucont_im_x;
        Ucont_ig_y = Ucont_im_y;

        Ucat_ig_x = Ucat_im_x;
        Ucat_ig_y = Ucat_im_y;
        
        pseudo_t = pseudo_t+1;     
        
        fprintf('subitr = %d, Momentum convergence error = %10e \n', pseudo_t, err); 
    end 
    
    if (err > 1e-1)
        fprintf('time = %d, Momentum solver does not converge e = %10e \n', time, err);
    end

end