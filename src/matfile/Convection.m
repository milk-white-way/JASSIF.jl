function [Conv_Flux_x Conv_Flux_y] = Convection(Ucont_x, Ucont_y, Ucat_x, Ucat_y,dx, dy)

% Fpx , Fpy are convective flux at the half node
%% Conv_Flux* are fluxes at the integer node
M = length(Ucat_x(:,1));
N = length(Ucat_x(1,:));
coef = 1/8;
% All the hafl node 
for i = 1:M-1
    for j = 1:N
        
         ucon = Ucont_x(i,j) /  2;
         
         up = ucon + abs(ucon);
         um = ucon - abs(ucon);
         
      if (i > 1 && i < M-1)
            Fpx1(i,j) =  um * ( coef * ( - Ucat_x(i+2,j) - 2 * Ucat_x(i+1,j) + 3 * Ucat_x(i,j)  ) + Ucat_x(i+1,j)) +  up * ( coef * ( - Ucat_x(i-1,j) - 2 * Ucat_x(i,j)   + 3 * Ucat_x(i+1,j)) + Ucat_x(i,j));                          
      else
          if (i == 1)
              Fpx1(i,j) =um * ( coef * ( - Ucat_x(i+2,j) - 2 * Ucat_x(i+1,j) + 3 * Ucat_x(i,j)  ) + Ucat_x(i+1,j)) +  up * ( coef * ( - Ucat_x(i,j)   - 2 * Ucat_x(i,j)   + 3 * Ucat_x(i+1,j)) + Ucat_x(i,j));                                                      
          else
                 if (i == M-1)
                    Fpx1(i,j) =  um * ( coef * ( - Ucat_x(i+1,j) - 2 * Ucat_x(i+1,j) + 3 * Ucat_x(i,j)  ) + Ucat_x(i+1,j)) + up * ( coef * ( - Ucat_x(i-1,j) - 2 * Ucat_x(i,j)   + 3 * Ucat_x(i+1,j)) + Ucat_x(i,j));                                             
                 end % End of i == M-1
          end % End of i == 1
      end % End of internal nodes
      
    end
end

% Zero-out the boundary
for j = 1:N    
    Fpx1(M,j) = 0;    
end

%% ----------------------------------
for i = 1:M
    for j = 1:N-1
        
      ucon = Ucont_y(i,j) /  2;
         
      up = ucon + abs(ucon);
      um = ucon - abs(ucon);
      
      if (j > 1 && j < N-1)
            Fpx2(i,j) = um * ( coef * ( - Ucat_x(i,j+2) - 2 * Ucat_x(i,j+1) + 3 * Ucat_x(i,j)  ) + Ucat_x(i,j+1)) +  up * ( coef * ( - Ucat_x(i,j-1) - 2 * Ucat_x(i,j)   + 3 * Ucat_x(i,j+1)) + Ucat_x(i,j));                          
      else
          if (j == 1)
              Fpx2(i,j) = um * ( coef * ( - Ucat_x(i,j+2) - 2 * Ucat_x(i,j+1) + 3 * Ucat_x(i,j)  ) + Ucat_x(i,j+1)) +  up * ( coef * ( - Ucat_x(i,j)   - 2 * Ucat_x(i,j)   + 3 * Ucat_x(i,j+1)) + Ucat_x(i,j));                                                     
          else
                 if (j == N-1)
                     Fpx2(i,j) =  um * ( coef * ( - Ucat_x(i,j+1) - 2 * Ucat_x(i,j+1) + 3 * Ucat_x(i,j)  ) + Ucat_x(i,j+1)) + up * ( coef * ( - Ucat_x(i,j-1) - 2 * Ucat_x(i,j)   + 3 * Ucat_x(i,j+1)) + Ucat_x(i,j));                                           
                 end % End of j == N-1
          end % End of j == 1
      end % End of internal nodes
      
    end
end
   
% Zero-out the boundary
for i = 1:M       
    Fpx2(i,N) = 0;        
end


%% ----------------------------------------- y - contribution--------------       
for i  = 1:M-1
    for j= 1:N
      
         ucon = Ucont_x(i,j) /  2;
         
         up = ucon + abs(ucon);
         um = ucon - abs(ucon);
      
      if (i > 1 && i < M-1)
            Fpy1(i,j) =  um * ( coef * ( - Ucat_y(i+2,j) - 2 * Ucat_y(i+1,j) + 3 * Ucat_y(i,j)  ) + Ucat_y(i+1,j)) + up * ( coef * ( - Ucat_y(i-1,j) - 2 * Ucat_y(i,j)   + 3 * Ucat_y(i+1,j)) + Ucat_y(i,j));                          
      else
          if (i == 1)
              Fpy1(i,j) =um * ( coef * ( - Ucat_y(i+2,j) - 2 * Ucat_y(i+1,j) + 3 * Ucat_y(i,j)  ) + Ucat_y(i+1,j)) + up * ( coef * ( - Ucat_y(i,j)   - 2 * Ucat_y(i,j)   + 3 * Ucat_y(i+1,j)) + Ucat_y(i,j));                                                      
          else
                 if (i == M-1)
                    Fpy1(i,j) =  um * ( coef * ( - Ucat_y(i+1,j) - 2 * Ucat_y(i+1,j) + 3 * Ucat_y(i,j)  ) + Ucat_y(i+1,j)) + up * ( coef * ( - Ucat_y(i-1,j) - 2 * Ucat_y(i,j)   + 3 * Ucat_y(i+1,j)) + Ucat_y(i,j));                                        
                 end % End of i == M-1
          end % End of i == 1
      end % End of internal nodes
    end
end

for j = 1:N
    Fpy1(M,j) = 0;    
end
%%
for i =1:M
    for j =1:N-1
        
         ucon = Ucont_y(i,j) /  2;
         
         up = ucon + abs(ucon);
         um = ucon - abs(ucon);
         
      
         % y - contribution       
      if (j > 1 && j < N-1)
            Fpy2(i,j) = um * ( coef * ( - Ucat_y(i,j+2) - 2 * Ucat_y(i,j+1) + 3 * Ucat_y(i,j)  ) + Ucat_y(i,j+1)) + up * ( coef * ( - Ucat_y(i,j-1) - 2 * Ucat_y(i,j)   + 3 * Ucat_y(i,j+1)) + Ucat_y(i,j));                          
      else
          if (j == 1)
              Fpy2(i,j) = um * ( coef * ( - Ucat_y(i,j+2) - 2 * Ucat_y(i,j+1) + 3 * Ucat_y(i,j)  ) + Ucat_y(i,j+1)) + up * ( coef * ( - Ucat_y(i,j)   - 2 * Ucat_y(i,j)   +  3 * Ucat_y(i,j+1)) + Ucat_y(i,j));                                                      
          else
                 if (j == N-1)
                     Fpy2(i,j) =  um * ( coef * ( - Ucat_y(i,j+1) - 2 * Ucat_y(i,j+1) + 3 * Ucat_y(i,j)  ) + Ucat_y(i,j+1)) + up * ( coef * ( - Ucat_y(i,j-1) - 2 * Ucat_y(i,j)   + 3 * Ucat_y(i,j+1)) + Ucat_y(i,j));                                        
                 end % End of j == N-1
          end % End of j == 1
      end % End of internal nodes     
      
    end % End of j
end % End of i

for i = 1:M
    Fpy2(i,N) = 0;    
end


% Take the derivatives first term
for i = 2 :M-1
    for j = 2:N-1
        Conv_Flux_x(i,j) = (Fpx1(i,j) - Fpx1(i-1,j)) / dx  +  (Fpx2(i,j) - Fpx2(i,j-1)) / dy;        
        Conv_Flux_y(i,j) = (Fpy1(i,j) - Fpy1(i-1,j)) / dx   +  (Fpy2(i,j) - Fpy2(i,j-1)) / dy; 
    end
end

% Boundary zeroout
for i = 1:M
    Conv_Flux_x(i,1) = 0;
    Conv_Flux_y(i,1) = 0;
    Conv_Flux_x(i,N) = 0;
    Conv_Flux_y(i,N) = 0;
end

for j = 1:N
    Conv_Flux_x(1,j) = 0;
    Conv_Flux_y(1,j) = 0;
    Conv_Flux_x(M,j) = 0;
    Conv_Flux_y(M,j) = 0;
end

         