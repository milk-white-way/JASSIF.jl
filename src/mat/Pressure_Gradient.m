function [P_Gradient_x, P_Gradient_y] = ...
    Pressure_Gradient(M2, N2, Pressure, ...
                      iphys, iphye, jphys, jphye, ...
                      dx, dy, DEBUG)

    if DEBUG
        fprintf('DEBUG: \tCalculation for Pressure Gradient\n');
    end

    P_Gradient_x = zeros(M2, N2);
    P_Gradient_y = zeros(M2, N2);

    for ii = iphys:iphye
        for jj = jphys:jphye

            %% ---------- x direction ----------
            %if ii == 1
            %    % Forward difference
            %    P_Gradient_x(jj, ii) = ( -3*Pressure(jj, ii) + 4*Pressure(jj, ii+1) - Pressure(jj, ii+2) ) / (2*dx);
            %elseif ii == M2 
            %    % Backward difference
            %    P_Gradient_x(jj, ii) = ( Pressure(jj, ii-2) - 4*Pressure(jj, ii-1) + 3*Pressure(jj, ii) ) / (2*dx);
            %else
            % Central difference
            P_Gradient_x(ii, jj) = ( Pressure(ii+1, jj) - Pressure(ii-1, jj)) / (2*dx);
            %end
            
            %% ---------- y direction ----------
            %if jj == 1
            %    % Forward difference
            %    P_Gradient_y(jj, ii) = ( -3*Pressure(jj, ii) + 4*Pressure(jj+1, ii) - Pressure(jj+2, ii) ) / (2*dy);
            %elseif jj == N2
            %    % Backward difference
            %    P_Gradient_y(jj, ii) = ( Pressure(jj-2, ii) - 4*Pressure(jj-1, ii) + 3*Pressure(jj, ii) ) / (2*dy);
            %else
            % Central difference
            P_Gradient_y(ii, jj) = ( Pressure(ii, jj+1) - Pressure(ii, jj-1)) / (2*dy);
            %end
        end
    end

end