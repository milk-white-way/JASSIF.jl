function [Conv_Flux_x, Conv_Flux_y] = ...
    Convective_Flux(M, N, M2, N2, M3, N3, Nghost, ...
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

    coef = 1/8;

    % All the hafl node fluxes
    %% x - contribution
    if DEBUG
        fprintf('Fpx1 at \n');
    end
    %% ----------------------------------
    for i = 1:M3
        for j = 1:N
            ucon = Ucont_x(i, j) /  2;
            
            up = ucon + abs(ucon); % == 0 if ucont_x(i, j) < 0
            um = ucon - abs(ucon); % == 0 if ucont_x(i, j) > 0

            if DEBUG
                fprintf('\t (%d, %d) has Ucont_x that ', i, j);
            end

            if (i == 1)
                Fpx1(i, j) = um * ( coef * ( - Ucat_x(i+G +2, j+G) - 2 * Ucat_x(i+G +1, j+G) + 3 * Ucat_x(i+G   , j+G) ) + Ucat_x(i+G +1, j+G) ) ...
                           + up * ( coef * ( - Ucat_x(i+G   , j+G) - 2 * Ucat_x(i+G   , j+G) + 3 * Ucat_x(i+G +1, j+G) ) + Ucat_x(i+G   , j+G) );
                
                if DEBUG
                    if up
                        fprintf('> 0 and used Ucat_x at (%d, %d), (%d, %d) \n', i+G+1, j+G, i+G, j+G);
                    elseif um
                        fprintf('< 0 and used Ucat_x at (%d, %d), (%d, %d), (%d, %d) \n', i+G+2, j+G, i+G+1, j+G, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            elseif (i == M3)
                Fpx1(i, j) =  um * ( coef * ( - Ucat_x(i+G +1, j+G) - 2 * Ucat_x(i+G +1, j+G) + 3 * Ucat_x(i+G   , j+G) ) + Ucat_x(i+G+1, j+G) ) ...
                            + up * ( coef * ( - Ucat_x(i+G -1, j+G) - 2 * Ucat_x(i+G   , j+G) + 3 * Ucat_x(i+G +1, j+G) ) + Ucat_x(i+G  , j+G) );
                        
                if DEBUG
                    if up 
                        fprintf('< 0 and used Ucat_x at (%d, %d), (%d, %d), (%d, %d) \n', i+G+1, j+G, i+G, j+G, i+G-1, j+G);
                    elseif um
                        fprintf('> 0 and used Ucat_x at (%d, %d), (%d, %d) \n', i+G+1, j+G, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            else
                Fpx1(i, j) =  um * ( coef * ( - Ucat_x(i+G +2, j+G) - 2 * Ucat_x(i+G +1, j+G) + 3 * Ucat_x(i+G   , j+G) ) + Ucat_x(i+G +1, j+G) ) ...
                            + up * ( coef * ( - Ucat_x(i+G -1, j+G) - 2 * Ucat_x(i+G   , j+G) + 3 * Ucat_x(i+G +1, j+G) ) + Ucat_x(i+G   , j+G) );

                if DEBUG
                    if up 
                        fprintf('< 0 and used Ucat_x at (%d, %d), (%d, %d), (%d, %d) \n', i+G+1, j+G, i+G, j+G, i+G-1, j+G);
                    elseif um
                        fprintf('< 0 and used Ucat_x at (%d, %d), (%d, %d), (%d, %d) \n', i+G+2, j+G, i+G+1, j+G, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            end
        end
    end

    if DEBUG
        fprintf('Fpx2 at \n');
    end
    %% ----------------------------------
    for i = 1:M
        for j = 1:N3
            
            ucon = Ucont_y(i, j) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            if DEBUG
                fprintf('\t (%d, %d) has Ucont_y that ', i, j);
            end
        
            if (j == 1)
                Fpx2(i, j) = um * ( coef * ( - Ucat_x(i+G, j+G +2) - 2 * Ucat_x(i+G, j+G +1) + 3 * Ucat_x(i+G, j+G  )  ) + Ucat_x(i+G, j+G +1) ) ...
                           + up * ( coef * ( - Ucat_x(i+G, j+G   ) - 2 * Ucat_x(i+G, j+G   ) + 3 * Ucat_x(i+G, j+G +1) ) + Ucat_x(i+G, j+G   ) );
                
                if DEBUG
                    if up
                        fprintf('> 0 and used Ucat_x at (%d, %d), (%d, %d) \n', i+G, j+G+1, i+G, j+G);
                    elseif um
                        fprintf('< 0 and used Ucat_x at (%d, %d), (%d, %d), (%d, %d) \n', i+G, j+G+2, i+G, j+G+1, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            elseif (j == N3)
                Fpx2(i, j) =  um * ( coef * ( - Ucat_x(i+G, j+G +1) - 2 * Ucat_x(i+G, j+G +1) + 3 * Ucat_x(i+G,   j+G)  ) + Ucat_x(i+G, j+G +1) ) ...
                            + up * ( coef * ( - Ucat_x(i+G, j+G -1) - 2 * Ucat_x(i+G, j+G   ) + 3 * Ucat_x(i+G, j+G +1) ) + Ucat_x(i+G, j+G   ) );
                        
                if DEBUG
                    if up 
                        fprintf('< 0 and used Ucat_x at (%d, %d), (%d, %d), (%d, %d) \n', i+G, j+G+1, i+G, j+G, i+G-1, j+G);
                    elseif um
                        fprintf('> 0 and used Ucat_x at (%d, %d), (%d, %d) \n', i+G, j+G+1, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            else
                Fpx2(i, j) = um * ( coef * ( - Ucat_x(i+G, j+G +2) - 2 * Ucat_x(i+G, j+G +1) + 3 * Ucat_x(i+G,   j+G)  ) + Ucat_x(i+G, j+G +1) ) ...
                           + up * ( coef * ( - Ucat_x(i+G, j+G -1) - 2 * Ucat_x(i+G, j+G   ) + 3 * Ucat_x(i+G, j+G +1) ) + Ucat_x(i+G, j+G   ) );

                if DEBUG
                    if up 
                        fprintf('< 0 and used Ucat_x at (%d, %d), (%d, %d), (%d, %d) \n', i+G, j+G+1, i+G, j+G, i+G-1, j+G);
                    elseif um
                        fprintf('> 0 and used Ucat_x at (%d, %d), (%d, %d), (%d, %d) \n', i+G, j+G+2, i+G, j+G+1, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            end
        end
    end
    
    %% y - contribution
    if DEBUG
        fprintf('Fpy1 at \n');
    end
    %% ----------------------------------
    for i = 1:M3
        for j = 1:N
            ucon = Ucont_x(i, j) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            if DEBUG
                fprintf('\t (%d, %d) has Ucont_x that ', i, j);
            end

            if (i == 1)
                Fpy1(i, j) = um * ( coef * ( - Ucat_y(i+G +2, j+G) - 2 * Ucat_y(i+G +1, j+G) + 3 * Ucat_y(i+G   , j+G) ) + Ucat_y(i+G +1, j+G) ) ...
                           + up * ( coef * ( - Ucat_y(i+G   , j+G) - 2 * Ucat_y(i+G   , j+G) + 3 * Ucat_y(i+G +1, j+G) ) + Ucat_y(i+G   , j+G) );
                
                if DEBUG
                    if up
                        fprintf('> 0 and used Ucat_y at (%d, %d), (%d, %d) \n', i+G+1, j+G, i+G, j+G);
                    elseif um
                        fprintf('< 0 and used Ucat_y at (%d, %d), (%d, %d), (%d, %d) \n', i+G+2, j+G, i+G+1, j+G, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            elseif (i == M3)
                Fyx1(i, j) =  um * ( coef * ( - Ucat_y(i+G +1, j+G) - 2 * Ucat_y(i+G +1, j+G) + 3 * Ucat_y(i+G   , j+G) ) + Ucat_y(i+G+1, j+G) ) ...
                            + up * ( coef * ( - Ucat_y(i+G -1, j+G) - 2 * Ucat_y(i+G   , j+G) + 3 * Ucat_y(i+G +1, j+G) ) + Ucat_y(i+G  , j+G) );
                        
                if DEBUG
                    if up 
                        fprintf('< 0 and used Ucat_y at (%d, %d), (%d, %d), (%d, %d) \n', i+G+1, j+G, i+G, j+G, i+G-1, j+G);
                    elseif um
                        fprintf('> 0 and used Ucat_y at (%d, %d), (%d, %d) \n', i+G+1, j+G, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            else
                Fyx1(i, j) =  um * ( coef * ( - Ucat_y(i+G +2, j+G) - 2 * Ucat_y(i+G +1, j+G) + 3 * Ucat_y(i+G   , j+G) ) + Ucat_y(i+G +1, j+G) ) ...
                            + up * ( coef * ( - Ucat_y(i+G -1, j+G) - 2 * Ucat_y(i+G   , j+G) + 3 * Ucat_y(i+G +1, j+G) ) + Ucat_y(i+G   , j+G) );

                if DEBUG
                    if up 
                        fprintf('< 0 and used Ucat_y at (%d, %d), (%d, %d), (%d, %d) \n', i+G+1, j+G, i+G, j+G, i+G-1, j+G);
                    elseif um
                        fprintf('< 0 and used Ucat_y at (%d, %d), (%d, %d), (%d, %d) \n', i+G+2, j+G, i+G+1, j+G, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            end
        end
    end

    if DEBUG
        fprintf('Fpy2 at \n');
    end
    %% ----------------------------------
    for i = 1:M
        for j = 1:N3
            
            ucon = Ucont_y(i, j) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            if DEBUG
                fprintf('\t (%d, %d) has Ucont_y that ', i, j);
            end
        
            if (j == 1)
                Fpy2(i, j) = um * ( coef * ( - Ucat_y(i+G, j+G +2) - 2 * Ucat_y(i+G, j+G +1) + 3 * Ucat_y(i+G, j+G  )  ) + Ucat_y(i+G, j+G +1) ) ...
                           + up * ( coef * ( - Ucat_y(i+G, j+G   ) - 2 * Ucat_y(i+G, j+G   ) + 3 * Ucat_y(i+G, j+G +1) ) + Ucat_y(i+G, j+G   ) );
                
                if DEBUG
                    if up
                        fprintf('> 0 and used Ucat_y at (%d, %d), (%d, %d) \n', i+G, j+G+1, i+G, j+G);
                    elseif um
                        fprintf('< 0 and used Ucat_y at (%d, %d), (%d, %d), (%d, %d) \n', i+G, j+G+2, i+G, j+G+1, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            elseif (j == N3)
                Fpy2(i, j) =  um * ( coef * ( - Ucat_y(i+G, j+G +1) - 2 * Ucat_y(i+G, j+G +1) + 3 * Ucat_y(i+G,   j+G)  ) + Ucat_y(i+G, j+G +1) ) ...
                            + up * ( coef * ( - Ucat_y(i+G, j+G -1) - 2 * Ucat_y(i+G, j+G   ) + 3 * Ucat_y(i+G, j+G +1) ) + Ucat_y(i+G, j+G   ) );
   
                if DEBUG
                    if up 
                        fprintf('< 0 and used Ucat_y at (%d, %d), (%d, %d), (%d, %d) \n', i+G, j+G+1, i+G, j+G, i+G-1, j+G);
                    elseif um
                        fprintf('> 0 and used Ucat_y at (%d, %d), (%d, %d) \n', i+G, j+G+1, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            else
                Fpy2(i, j) = um * ( coef * ( - Ucat_y(i+G, j+G +2) - 2 * Ucat_y(i+G, j+G +1) + 3 * Ucat_y(i+G,   j+G)  ) + Ucat_y(i+G, j+G +1) ) ...
                           + up * ( coef * ( - Ucat_y(i+G, j+G -1) - 2 * Ucat_y(i+G, j+G   ) + 3 * Ucat_y(i+G, j+G +1) ) + Ucat_y(i+G, j+G   ) );

                if DEBUG
                    if up 
                        fprintf('< 0 and used Ucat_y at (%d, %d), (%d, %d), (%d, %d) \n', i+G, j+G+1, i+G, j+G, i+G-1, j+G);
                    elseif um
                        fprintf('> 0 and used Ucat_y at (%d, %d), (%d, %d), (%d, %d) \n', i+G, j+G+2, i+G, j+G+1, i+G, j+G);
                    else
                        fprintf('= 0 and thus, = 0 \n');
                    end
                end
            end
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