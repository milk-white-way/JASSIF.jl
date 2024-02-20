function [P_Gradient_x, P_Gradient_y] = ...
    Pressure_Gradient(FluxSumOld, Pressure, M2, N2, dx, dy)

    P_Gradient_x = FluxSumOld.P_Gradient.Flux_x;
    P_Gradient_y = FluxSumOld.P_Gradient.Flux_y;

    for ii = 1:M2
        for jj = 1:N2

            %% ---------- x direction ----------
            if ii == 1
                % Forward difference
                P_Gradient_x(jj, ii) = ( -3*Pressure(jj, ii) + 4*Pressure(jj, ii+1) - Pressure(jj, ii+2) ) / (2*dx);
            elseif ii == M2 
                % Backward difference
                P_Gradient_x(jj, ii) = ( Pressure(jj, ii-2) - 4*Pressure(jj, ii-1) + 3*Pressure(jj, ii) ) / (2*dx);
            else
                % Central difference
                P_Gradient_x(jj, ii) = ( Pressure(jj, ii+1) - Pressure(jj, ii-1)) / (2*dx);
            end
            
            %% ---------- y direction ----------
            if jj == 1
                % Forward difference
                P_Gradient_y(jj, ii) = ( -3*Pressure(jj, ii) + 4*Pressure(jj+1, ii) - Pressure(jj+2, ii) ) / (2*dy);
            elseif jj == N2
                % Backward difference
                P_Gradient_y(jj, ii) = ( Pressure(jj-2, ii) - 4*Pressure(jj-1, ii) + 3*Pressure(jj, ii) ) / (2*dy);
            else
                % Central difference
                P_Gradient_y(jj, ii) = ( Pressure(jj+1, ii) - Pressure(jj-1, ii)) / (2*dy);
            end
        end
    end

end