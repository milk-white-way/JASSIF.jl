function [Viscous_x, Viscous_y] = ...
    Viscous_Flux(M2, N2, Ucat_x, Ucat_y, ...
                 iphys, iphye, jphys, jphye, ...
                 dx, dy, Re, DEBUG) % Re is Reynolds number

    if DEBUG
        fprintf('DEBUG: \tCalculation of Fluxes for Viscous Terms\n');
    end

    Viscous_x = zeros(M2, N2);
    Viscous_y = zeros(M2, N2);

    %% Inner (Physical domain)
    for ii = iphys:iphye
        for jj = jphys:jphye

            %% ------------- x direction -------------------
            % Central differencing
            up = Ucat_x(ii, jj);
            % E-W
            ue = Ucat_x(ii+1, jj);
            uw = Ucat_x(ii-1, jj);
            % N-S
            un = Ucat_x(ii, jj+1);
            us = Ucat_x(ii, jj-1);
                
            Viscous_x(ii, jj) = ( (uw - 2*up + ue)/(dx^2) + (us - 2*up + un)/(dy^2) );

            %% ------------- y direction --------------------
            % Central differencing
            vp = Ucat_y(ii, jj);
            % E-W
            ve = Ucat_y(ii+1, jj);
            vw = Ucat_y(ii-1, jj);
            % N-S
            vn = Ucat_y(ii, jj+1);
            vs = Ucat_y(ii, jj-1);
                
            Viscous_y(ii, jj) = ( (vw - 2*vp + ve)/(dx^2) + (vs - 2*vp + vn)/(dy^2) );

        end
    end

    %% Outer: does not exists in the context of periodic boundary conditions
    %{
    % South (jj = 1) and North (jj = N2)
    for jj = [1, N2]
        for ii = 2:M2-1 % Minus the corners
            if jj == 1 % Forwards differencing
                if DEBUG 
                    fprintf('DEBUG: \tSouth bound at (%d, %d)\n', ii, jj);
                end
                %% ------------- x direction -------------------
                up = Ucat_x(jj, ii);
                % E-W
                ue = Ucat_x(jj, ii+1);
                uw = Ucat_x(jj, ii-1);
                % N-S
                udy = Ucat_x(jj+1, ii);
                uddy = Ucat_x(jj+2, ii);
                udddy = Ucat_x(jj+3, ii);

                Viscous_x(jj, ii) = ( (uw - 2*up + ue)/(dx^2) + (2*up - 5*udy + 4*uddy - udddy)/(dy^3) );

                %% ------------- y direction --------------------
                vp = Ucat_y(jj, ii);
                % E-W
                ve = Ucat_y(jj, ii+1);
                vw = Ucat_y(jj, ii-1);
                % N-S
                vdy = Ucat_y(jj+1, ii);
                vddy = Ucat_y(jj+2, ii);
                vdddy = Ucat_y(jj+3, ii);
                
                Viscous_y(jj, ii) = ( (vw - 2*vp + ve)/(dx^2) + (2*vp - 5*vdy + 4*vddy - vdddy)/(dy^3) );

            elseif jj == N2 % Backwards differencing
                if DEBUG 
                    fprintf('DEBUG: \tNorth bound at (%d, %d)\n', ii, jj);
                end
                %% ------------- x direction -------------------
                up = Ucat_x(jj, ii);
                % E-W
                ue = Ucat_x(jj, ii+1);
                uw = Ucat_x(jj, ii-1);
                % N-S
                udy = Ucat_x(jj-1, ii);
                uddy = Ucat_x(jj-2, ii);
                udddy = Ucat_x(jj-3, ii);

                Viscous_x(jj, ii) = ( (uw - 2*up + ue)/(dx^2) + (2*up - 5*udy + 4*uddy - udddy)/(dy^3) );

                %% ------------- y direction --------------------
                vp = Ucat_y(jj, ii);
                % E-W
                ve = Ucat_y(jj, ii+1);
                vw = Ucat_y(jj, ii-1);
                % N-S
                vdy = Ucat_y(jj-1, ii);
                vddy = Ucat_y(jj-2, ii);
                vdddy = Ucat_y(jj-3, ii);
                
                Viscous_y(jj, ii) = ( (vw - 2*vp + ve)/(dx^2) + (2*vp - 5*vdy + 4*vddy - vdddy)/(dy^3) );

            end
        end
    end
    
    % West (ii = 1) and East (ii = M2)
    for jj = 2:N2-1
        for ii = [1, M2] % Minus the corners
            if ii == 1 % Forwards differencing
                if DEBUG 
                    fprintf('DEBUG: \tWest bound at (%d, %d)\n', ii, jj);
                end
                %% ------------- x direction -------------------
                % Forwards differencing
                up = Ucat_x(jj, ii);
                % E-W
                udx = Ucat_x(jj, ii+1);
                uddx = Ucat_x(jj, ii+2);
                udddx = Ucat_x(jj, ii+3);
                % N-S
                un = Ucat_x(jj+1, ii);
                us = Ucat_x(jj-1, ii);

                Viscous_x(jj, ii) = ( (2*up - 5*udx + 4*uddx - udddx)/(dx^3) + (us - 2*up + un)/(dy^2) );

                %% ------------- y direction --------------------
                vp = Ucat_y(jj, ii);
                % E-W
                vdx = Ucat_y(jj, ii+1);
                vddx = Ucat_y(jj, ii+2);
                vdddx = Ucat_y(jj, ii+3);
                % N-S
                vn = Ucat_y(jj+1, ii);
                vs = Ucat_y(jj-1, ii);

                Viscous_y(jj, ii) = ( (2*vp - 5*vdx + 4*vddx - vdddx)/(dx^3) + (vs - 2*vp + vn)/(dy^2) );

            elseif ii == M2 % Backwards differencing
                if DEBUG 
                    fprintf('DEBUG: \tEast bound at (%d, %d)\n', ii, jj);
                end
                %% ------------- x direction -------------------
                % Forwards differencing
                up = Ucat_x(jj, ii);
                % E-W
                udx = Ucat_x(jj, ii-1);
                uddx = Ucat_x(jj, ii-2);
                udddx = Ucat_x(jj, ii-3);
                % N-S
                un = Ucat_x(jj+1, ii);
                us = Ucat_x(jj-1, ii);

                Viscous_x(jj, ii) = ( (2*up - 5*udx + 4*uddx - udddx)/(dx^3) + (us - 2*up + un)/(dy^2) );

                %% ------------- y direction --------------------
                vp = Ucat_y(jj, ii);
                % E-W
                vdx = Ucat_y(jj, ii-1);
                vddx = Ucat_y(jj, ii-2);
                vdddx = Ucat_y(jj, ii-3);
                % N-S
                vn = Ucat_y(jj+1, ii);
                vs = Ucat_y(jj-1, ii);

                Viscous_y(jj, ii) = ( (2*vp - 5*vdx + 4*vddx - vdddx)/(dx^3) + (vs - 2*vp + vn)/(dy^2) );

            end
        end
    end

    % Corners
    % |
    % |
    % o------
    Viscous_x(1, 1) = ( (2*Ucat_x(1, 1) - 5*Ucat_x(1, 2) + 4*Ucat_x(1, 3) - Ucat_x(1, 4))/(dx^3) + (2*Ucat_x(1, 1) - 5*Ucat_x(2, 1) + 4*Ucat_x(3, 1) - Ucat_x(4, 1))/(dy^3) );
    Viscous_y(1, 1) = ( (2*Ucat_y(1, 1) - 5*Ucat_y(1, 2) + 4*Ucat_y(1, 3) - Ucat_y(1, 4))/(dx^3) + (2*Ucat_y(1, 1) - 5*Ucat_y(2, 1) + 4*Ucat_y(3, 1) - Ucat_y(4, 1))/(dy^3) );

    %       |
    %       |
    % ------o
    Viscous_x(1, M2) = ( (2*Ucat_x(1, M2) - 5*Ucat_x(1, M2-1) + 4*Ucat_x(1, M2-2) - Ucat_x(1, M2-3))/(dx^3) + (2*Ucat_x(1, M2) - 5*Ucat_x(2, M2) + 4*Ucat_x(3, M2) - Ucat_x(4, M2))/(dy^3) );
    Viscous_y(1, M2) = ( (2*Ucat_y(1, M2) - 5*Ucat_y(1, M2-1) + 4*Ucat_y(1, M2-2) - Ucat_y(1, M2-3))/(dx^3) + (2*Ucat_y(1, M2) - 5*Ucat_y(2, M2) + 4*Ucat_y(3, M2) - Ucat_y(4, M2))/(dy^3) );

    % o------
    % |
    % |
    Viscous_x(N2, 1) = ( (2*Ucat_x(N2, 1) - 5*Ucat_x(N2, 2) + 4*Ucat_x(N2, 3) - Ucat_x(N2, 4))/(dx^3) + (2*Ucat_x(N2, 1) - 5*Ucat_x(N2-1, 1) + 4*Ucat_x(N2-2, 1) - Ucat_x(N2-3, 1))/(dy^3) );
    Viscous_y(N2, 1) = ( (2*Ucat_y(N2, 1) - 5*Ucat_y(N2, 2) + 4*Ucat_y(N2, 3) - Ucat_y(N2, 4))/(dx^3) + (2*Ucat_y(N2, 1) - 5*Ucat_y(N2-1, 1) + 4*Ucat_y(N2-2, 1) - Ucat_y(N2-3, 1))/(dy^3) );

    % ------o
    %       |
    %       |
    Viscous_x(N2, M2) = ( (2*Ucat_x(N2, M2) - 5*Ucat_x(N2, M2-1) + 4*Ucat_x(N2, M2-2) - Ucat_x(N2, M2-3))/(dx^3) + (2*Ucat_x(N2, M2) - 5*Ucat_x(N2-1, M2) + 4*Ucat_x(N2-2, M2) - Ucat_x(N2-3, M2))/(dy^3) );
    Viscous_y(N2, M2) = ( (2*Ucat_y(N2, M2) - 5*Ucat_y(N2, M2-1) + 4*Ucat_y(N2, M2-2) - Ucat_y(N2, M2-3))/(dx^3) + (2*Ucat_y(N2, M2) - 5*Ucat_y(N2-1, M2) + 4*Ucat_y(N2-2, M2) - Ucat_y(N2-3, M2))/(dy^3) );

    %}

    Viscous_x = Viscous_x/Re;
    Viscous_y = Viscous_y/Re;

end