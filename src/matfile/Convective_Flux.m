function [Conv_Flux_x, Conv_Flux_y] = ...
    Convective_Flux(FluxSumOld, ...
                        Ucont_x, Ucont_y, Ucat_x, Ucat_y, ...
                        M, N, Nghost, dx, dy, ...
                        iphys, iphye, jphys, jphye, DEBUG)

    Conv_Flux_x = FluxSumOld.Convective.Flux_x;
    Conv_Flux_y = FluxSumOld.Convective.Flux_y;

    M3 = M+1;
    N3 = N+1;

    G = Nghost;

    coef = 1/8;

    for ii = iphys:iphye
        for jj = jphys:jphye
            Conv_Flux_y(jj, ii) = 0;
            Conv_Flux_x(jj, ii) = 0;
        end
    end

    % All the hafl node fluxes
    %% x - contribution
    if DEBUG
        fprintf('Fpx1 at \n');
    end
    %% ----------------------------------
    for i = 1:M3
        for j = 1:N
            ucon = Ucont_x(j, i) /  2;
            
            up = ucon + abs(ucon); % == 0 if ucont_x(j, i) < 0
            um = ucon - abs(ucon); % == 0 if ucont_x(j, i) > 0

            if DEBUG
                fprintf('\t (%d, %d) has Ucont_x that ', i, j);
            end

            if (i == 1)
                Fpx1(j, i) = um * ( coef * ( - Ucat_x(j+G, i+G+2) - 2 * Ucat_x(j+G, i+G+1) + 3 * Ucat_x(j+G, i+G)  ) + Ucat_x(j+G, i+G+1) ) ...
                           + up * ( coef * ( - Ucat_x(j+G, i+G)   - 2 * Ucat_x(j+G, i+G)   + 3 * Ucat_x(j+G, i+G+1)) + Ucat_x(j+G, i+G)   );
                
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
                Fpx1(j, i) =  um * ( coef * ( - Ucat_x(j+G, i+G+1) - 2 * Ucat_x(j+G, i+G+1) + 3 * Ucat_x(j+G, i+G)  ) + Ucat_x(j+G, i+G+1) ) ...
                            + up * ( coef * ( - Ucat_x(j+G, i+G-1) - 2 * Ucat_x(j+G, i+G)   + 3 * Ucat_x(j+G, i+G+1)) + Ucat_x(j+G, i+G)   );
                        
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
                Fpx1(j, i) =  um * ( coef * ( - Ucat_x(j+G, i+G+2) - 2 * Ucat_x(j+G, i+G+1) + 3 * Ucat_x(j+G, i+G)  ) + Ucat_x(j+G, i+G+1) ) ...
                            + up * ( coef * ( - Ucat_x(j+G, i+G-1) - 2 * Ucat_x(j+G, i+G)   + 3 * Ucat_x(j+G, i+G+1)) + Ucat_x(j+G, i+G)   );

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
            
            ucon = Ucont_y(j, i) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            if DEBUG
                fprintf('\t (%d, %d) has Ucont_y that ', i, j);
            end
        
            if (j == 1)
                Fpx2(j, i) = um * ( coef * ( - Ucat_x(j+G+2, i+G) - 2 * Ucat_x(j+G+1, i+G) + 3 * Ucat_x(j+G,   i+G) ) + Ucat_x(j+G+1, i+G) ) ...
                           + up * ( coef * ( - Ucat_x(j+G,   i+G) - 2 * Ucat_x(j+G,   i+G) + 3 * Ucat_x(j+G+1, i+G))  + Ucat_x(j+G,   i+G) );
                
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
                Fpx2(j, i) =  um * ( coef * ( - Ucat_x(j+G+1, i+G) - 2 * Ucat_x(j+G+1, i+G) + 3 * Ucat_x(j+G,   i+G) ) + Ucat_x(j+G+1, i+G) ) ...
                            + up * ( coef * ( - Ucat_x(j+G-1, i+G) - 2 * Ucat_x(j+G,   i+G) + 3 * Ucat_x(j+G+1, i+G))  + Ucat_x(j+G,   i+G) );
                        
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
                Fpx2(j, i) = um * ( coef * ( - Ucat_x(j+G+2, i+G) - 2 * Ucat_x(j+G+1, i+G) + 3 * Ucat_x(j+G,   i+G) ) + Ucat_x(j+G+1, i+G) ) ...
                           + up * ( coef * ( - Ucat_x(j+G-1, i+G) - 2 * Ucat_x(j+G,   i+G) + 3 * Ucat_x(j+G+1, i+G))  + Ucat_x(j+G,   i+G) );

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
            ucon = Ucont_x(j, i) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            if DEBUG
                fprintf('\t (%d, %d) has Ucont_x that ', i, j);
            end

            if (i == 1)
                Fpy1(j, i) = um * ( coef * ( - Ucat_y(j+G, i+G+2) - 2 * Ucat_y(j+G, i+G+1) + 3 * Ucat_y(j+G, i+G)  ) + Ucat_y(j+G, i+G+1) ) ...
                           + up * ( coef * ( - Ucat_y(j+G, i+G)   - 2 * Ucat_y(j+G, i+G)   + 3 * Ucat_y(j+G, i+G+1)) + Ucat_y(j+G, i+G)   );
                
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
                Fpy1(j, i) =  um * ( coef * ( - Ucat_y(j+G, i+G+1) - 2 * Ucat_y(j+G, i+G+1) + 3 * Ucat_y(j+G, i+G)  ) + Ucat_y(j+G, i+G+1) ) ...
                            + up * ( coef * ( - Ucat_y(j+G, i+G-1) - 2 * Ucat_y(j+G, i+G)   + 3 * Ucat_y(j+G, i+G+1)) + Ucat_y(j+G, i+G)   );
                        
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
                Fpy1(j, i) =  um * ( coef * ( - Ucat_y(j+G, i+G+2) - 2 * Ucat_y(j+G, i+G+1) + 3 * Ucat_y(j+G, i+G)  ) + Ucat_y(j+G, i+G+1) ) ...
                            + up * ( coef * ( - Ucat_y(j+G, i+G-1) - 2 * Ucat_y(j+G, i+G)   + 3 * Ucat_y(j+G, i+G+1)) + Ucat_y(j+G, i+G)   );

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
            
            ucon = Ucont_y(j, i) /  2;
            
            up = ucon + abs(ucon);
            um = ucon - abs(ucon);

            if DEBUG
                fprintf('\t (%d, %d) has Ucont_y that ', i, j);
            end
        
            if (j == 1)
                Fpy2(j, i) = um * ( coef * ( - Ucat_y(j+G+2, i+G) - 2 * Ucat_y(j+G+1, i+G) + 3 * Ucat_y(j+G,   i+G) ) + Ucat_y(j+G+1, i+G) ) ...
                           + up * ( coef * ( - Ucat_y(j+G,   i+G) - 2 * Ucat_y(j+G,   i+G) + 3 * Ucat_y(j+G+1, i+G))  + Ucat_y(j+G,   i+G) );
                
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
                Fpy2(j, i) =  um * ( coef * ( - Ucat_y(j+G+1, i+G) - 2 * Ucat_y(j+G+1, i+G) + 3 * Ucat_y(j+G,   i+G) ) + Ucat_y(j+G+1, i+G) ) ...
                            + up * ( coef * ( - Ucat_y(j+G-1, i+G) - 2 * Ucat_y(j+G,   i+G) + 3 * Ucat_y(j+G+1, i+G))  + Ucat_y(j+G,   i+G) );
   
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
                Fpy2(j, i) = um * ( coef * ( - Ucat_y(j+G+2, i+G) - 2 * Ucat_y(j+G+1, i+G) + 3 * Ucat_y(j+G,   i+G) ) + Ucat_y(j+G+1, i+G) ) ...
                           + up * ( coef * ( - Ucat_y(j+G-1, i+G) - 2 * Ucat_y(j+G,   i+G) + 3 * Ucat_y(j+G+1, i+G))  + Ucat_y(j+G,   i+G) );

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
    for i = iphys:iphye
        for j = jphys:jphye
            Conv_Flux_x(j, i) = (Fpx1(j-G, i-G+1) - Fpx1(j-G, i-G)) / dx  +  (Fpx2(j-G+1, i-G) - Fpx2(j-G, i-G)) / dy; 
            Conv_Flux_y(j, i) = (Fpy1(j-G, i-G+1) - Fpy1(j-G, i-G)) / dx  +  (Fpy2(j-G+1, i-G) - Fpy2(j-G, i-G)) / dy;
        end
    end

end