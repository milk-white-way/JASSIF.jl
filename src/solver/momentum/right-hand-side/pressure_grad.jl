module PressureGradient
   export rhs_pressure_gradient_calc

   function rhs_pressure_gradient_calc(CSolution)
      let M = deepcopy(CSolution.m2)
          N = deepcopy(CSolution.n2)
          dx = deepcopy(CSolution.dx)
          dy = deepcopy(CSolution.dy)
          pres = deepcopy(CSolution.pres)

          pres_grad_x = Matrix{Float64}(undef, M, N)
          pres_grad_y = Matrix{Float64}(undef, M, N)

         for ii = 1:M
            for jj = 1:N
               if ( ii == 1 ) 
                  pres_grad_x[ii, jj] = 0
               elseif ( ii == 2 )
                  pres_grad_x[ii, jj] = ( pres[ii+1, jj] - pres[ii, jj] )/dx
               elseif ( ii == M-1 )
                  pres_grad_x[ii, jj] = ( pres[ii, jj] - pres[ii-1, jj] )/dx
               elseif ( ii == M )
                  pres_grad_x[ii, jj] = 0
               else
                  pres_grad_x[ii, jj] = ( pres[ii+1, jj] - pres[ii-1, jj] )/( 2*dx )
               end

               if ( jj == 1 ) 
                  pres_grad_y[ii, jj] = 0
               elseif ( jj == 2 )
                  pres_grad_y[ii, jj] = ( pres[ii, jj+1] - pres[ii, jj] )/dy
               elseif ( jj == N-1 )
                  pres_grad_y[ii, jj] = ( pres[ii, jj] - pres[ii, jj-1] )/dy
               elseif ( jj == N )
                  pres_grad_y[ii, jj] = 0
               else
                  pres_grad_y[ii, jj] = ( pres[ii, jj+1] - pres[ii, jj-1] )/( 2*dy )
               end
            end
         end

         return pres_grad_x, pres_grad_y
      end
   end
end