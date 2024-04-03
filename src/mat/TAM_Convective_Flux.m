function [Conv_Flux_x, Conv_Flux_y] = ...
    TAM_Convective_Flux(M, N, M2, N2, M3, N3, Nghost, ...
                    Ucont_x, Ucont_y, Ucat_x, Ucat_y, ...
                    iphys, iphye, jphys, jphye, ...
                    dx, dy, DEBUG)

    if DEBUG
        fprintf('DEBUG: \tCalculation of Fluxes for Convective Terms\n');
    end

    Fpx1 = zeros(M3, N);
    Fpx2 = zeros(M, N3);

    Fpy1 = zeros(M3, N);
    Fpy2 = zeros(M, N3);

    Conv_Flux_x = zeros(M2, N2);
    Conv_Flux_y = zeros(M2, N2);

    G = Nghost;

    alp = 1/2;
    bet = 1/8;

    % All the hafl node fluxes
    %% x - contribution
    %% ----------------------------------
    for ii = 1:M3
        for jj = 1:N
            ucon = Ucont_x(ii, jj) /  2;
            
            up = ucon + abs(ucon); % == 0 if ucont_x(ii, jj) < 0
            um = ucon - abs(ucon); % == 0 if ucont_x(ii, jj) > 0

            % west face for QUICK scheme
            Ucat_x_WW = Ucat_x(ii+G-2, jj+G);
            Ucat_x_W  = Ucat_x(ii+G-1, jj+G);
            Ucat_x_P  = Ucat_x(ii+G  , jj+G);
            Ucat_x_E  = Ucat_x(ii+G+1, jj+G);

            Fpx1(ii, jj) = up * ( alp*( Ucat_x_P + Ucat_x_W ) - bet*( Ucat_x_WW - 2*Ucat_x_W + Ucat_x_P ) ) ...
                         + um * ( alp*( Ucat_x_P + Ucat_x_W ) - bet*( Ucat_x_W  - 2*Ucat_x_P + Ucat_x_E ) );

        end
    end

    %% ----------------------------------
    for ii = 1:M
        for jj = 1:N3
            ucon = Ucont_y(ii,jj) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            % south face for QUICK scheme
            Ucat_x_SS = Ucat_x(ii+G, jj+G-2);
            Ucat_x_S  = Ucat_x(ii+G, jj+G-1);
            Ucat_x_P  = Ucat_x(ii+G, jj+G  );
            Ucat_x_N  = Ucat_x(ii+G, jj+G+1);

            Fpx2(ii, jj) = um * ( alp*( Ucat_x_P + Ucat_x_S ) - bet*( Ucat_x_SS - 2*Ucat_x_S + Ucat_x_P ) ) ...
                         + up * ( alp*( Ucat_x_P + Ucat_x_S ) - bet*( Ucat_x_S  - 2*Ucat_x_P + Ucat_x_N ) );

        end
    end
    
    %% y - contribution
    %% ----------------------------------
    for ii = 1:M3
        for jj = 1:N
            ucon = Ucont_x(ii, jj) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            % west face for QUICK scheme
            Ucat_y_WW = Ucat_y(ii+G-2, jj+G);
            Ucat_y_W  = Ucat_y(ii+G-1, jj+G);
            Ucat_y_P  = Ucat_y(ii+G  , jj+G);
            Ucat_y_E  = Ucat_y(ii+G+1, jj+G);

            Fpy1(ii, jj) = up * ( alp*( Ucat_y_P + Ucat_y_W ) - bet*( Ucat_y_WW - 2*Ucat_y_W + Ucat_y_P ) ) ...
                         + um * ( alp*( Ucat_y_P + Ucat_y_W ) - bet*( Ucat_y_W  - 2*Ucat_y_P + Ucat_y_E ) );

        end
    end

    %% ----------------------------------
    for ii = 1:M
        for jj = 1:N3
            
            ucon = Ucont_y(ii, jj) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            % south face for QUICK scheme
            Ucat_y_SS = Ucat_y(ii+G, jj+G-2);
            Ucat_y_S  = Ucat_y(ii+G, jj+G-1);
            Ucat_y_P  = Ucat_y(ii+G, jj+G  );
            Ucat_y_N  = Ucat_y(ii+G, jj+G+1);

            Fpy2(ii, jj) = um * ( alp*( Ucat_y_P + Ucat_y_S ) - bet*( Ucat_y_SS - 2*Ucat_y_S + Ucat_y_P ) ) ...
                         + up * ( alp*( Ucat_y_P + Ucat_y_S ) - bet*( Ucat_y_S  - 2*Ucat_y_P + Ucat_y_N ) );

        end
    end

    % Take the derivatives first term
    for ii = iphys:iphye
        for jj = jphys:jphye
            Conv_Flux_x(ii, jj) = ( Fpx1(ii-G +1, jj-G) - Fpx1(ii-G, jj-G) )/dx  +  ( Fpx2(ii-G, jj-G +1) - Fpx2(ii-G, jj-G) )/dy; 
            Conv_Flux_y(ii, jj) = ( Fpy1(ii-G +1, jj-G) - Fpy1(ii-G, jj-G) )/dx  +  ( Fpy2(ii-G, jj-G +1) - Fpy2(ii-G, jj-G) )/dy;
        end
    end

end