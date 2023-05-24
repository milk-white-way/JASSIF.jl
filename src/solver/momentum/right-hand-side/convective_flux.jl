module ConvectiveFlux
   export rhs_convective_flux_calc

   function rhs_convective_flux_calc(CSolution, coeff)
      let M = deepcopy(CSolution.m2)
          N = deepcopy(CSolution.n2)
          dx = deepcopy(CSolution.dx)
          dy = deepcopy(CSolution.dy)
          ucat_x = deepcopy(CSolution.ucat_x)
          ucat_y = deepcopy(CSolution.ucat_y)
          ucur_x = deepcopy(CSolution.ucur_x)
          ucur_y = deepcopy(CSolution.ucur_y)

          fpx1 = Matrix{Float64}(undef, M, N)
          fpy1 = Matrix{Float64}(undef, M, N)
          fpx2 = Matrix{Float64}(undef, M, N)
          fpy2 = Matrix{Float64}(undef, M, N)

          conv_flux_x = Matrix{Float64}(undef, M, N)
          conv_flux_y = Matrix{Float64}(undef, M, N)

         for ii = 1:M-1
            for jj = 1:N
               un = 0.5*ucur_x[ii, jj]
               um = un - abs(un)
               up = un + abs(un)

               vn = 0.5*ucur_y[ii, jj]
               vm = vn - abs(vn)
               vp = vn + abs(vn)

               if ( ii == 1 )
                  fpx1[ii,  jj] = um*( coeff*( - ucat_x[ii+2, jj] - 2*ucat_x[ii+1, jj] + 3*ucat_x[ii, jj] ) + ucat_x[ii+1, jj] ) + up*( coeff*( - ucat_x[ii, jj] - 2*ucat_x[ii, jj] + 3*ucat_x[ii+1, jj] ) + ucat_x[ii, jj] )
                  
                  fpy1[ii,  jj] = vm*( coeff*( - ucat_y[ii+2, jj] - 2*ucat_y[ii+1, jj] + 3*ucat_y[ii, jj] ) + ucat_y[ii+1, jj] ) + vp*( coeff*( - ucat_y[ii, jj] - 2*ucat_y[ii, jj] + 3*ucat_y[ii+1, jj] ) + ucat_y[ii, jj] )
               elseif ( ii == M-1 )
                  fpx1[ii,  jj] = um*( coeff*( - ucat_x[ii+1, jj] - 2*ucat_x[ii+1, jj] + 3*ucat_x[ii, jj] ) + ucat_x[ii+1, jj] ) + up*( coeff*( - ucat_x[ii-1, jj] - 2*ucat_x[ii, jj] + 3*ucat_x[ii+1, jj] ) + ucat_x[ii, jj] )
                  
                  fpy1[ii,  jj] = vm*( coeff*( - ucat_y[ii+1, jj] - 2*ucat_y[ii+1, jj] + 3*ucat_y[ii, jj] ) + ucat_y[ii+1, jj] ) + vp*( coeff*( - ucat_y[ii-1, jj] - 2*ucat_y[ii, jj] + 3*ucat_y[ii+1, jj] ) + ucat_y[ii, jj] )
               else
                  fpx1[ii,  jj] = um*( coeff*( - ucat_x[ii+2, jj] - 2*ucat_x[ii+1, jj] + 3*ucat_x[ii, jj] ) + ucat_x[ii+1, jj] ) + up*( coeff*( - ucat_x[ii-1, jj] - 2*ucat_x[ii, jj] + 3*ucat_x[ii+1, jj] ) + ucat_x[ii, jj] )
                  
                  fpy1[ii,  jj] = vm*( coeff*( - ucat_y[ii+2, jj] - 2*ucat_y[ii+1, jj] + 3*ucat_y[ii, jj] ) + ucat_y[ii+1, jj] ) + vp*( coeff*( - ucat_y[ii-1, jj] - 2*ucat_y[ii, jj] + 3*ucat_y[ii+1, jj] ) + ucat_y[ii, jj] )
               end
            end
         end

         for jj = 1:N
            fpx1[M, jj] = 0
            fpy1[M, jj] = 0
         end

         for ii = 1:M
            for jj = 1:N-1
               un = 0.5*ucur_x[ii, jj]
               um = un - abs(un)
               up = un + abs(un)

               vn = 0.5*ucur_y[ii, jj]
               vm = vn - abs(vn)
               vp = vn + abs(vn)

               if ( jj == 1 )
                  fpx2[ii,  jj] = um*( coeff*( - ucat_x[ii, jj+2] - 2*ucat_x[ii, jj+1] + 3*ucat_x[ii, jj] ) + ucat_x[ii, jj+1] ) + up*( coeff*( - ucat_x[ii, jj] - 2*ucat_x[ii, jj] + 3*ucat_x[ii, jj+1] ) + ucat_x[ii, jj] )
                  
                  fpy2[ii,  jj] = vm*( coeff*( - ucat_y[ii, jj+2] - 2*ucat_y[ii, jj+1] + 3*ucat_y[ii, jj] ) + ucat_y[ii, jj+1] ) + vp*( coeff*( - ucat_y[ii, jj] - 2*ucat_y[ii, jj] + 3*ucat_y[ii, jj+1] ) + ucat_y[ii, jj] )
               elseif ( jj == N-1 )
                  fpx2[ii,  jj] = um*( coeff*( - ucat_x[ii, jj+1] - 2*ucat_x[ii, jj+1] + 3*ucat_x[ii, jj] ) + ucat_x[ii, jj+1] ) + up*( coeff*( - ucat_x[ii, jj-1] - 2*ucat_x[ii, jj] + 3*ucat_x[ii, jj+1] ) + ucat_x[ii, jj] )
                  
                  fpy2[ii,  jj] = vm*( coeff*( - ucat_y[ii, jj+1] - 2*ucat_y[ii, jj+1] + 3*ucat_y[ii, jj] ) + ucat_y[ii, jj+1] ) + vp*( coeff*( - ucat_y[ii, jj-1] - 2*ucat_y[ii, jj] + 3*ucat_y[ii, jj+1] ) + ucat_y[ii, jj] )
               else
                  fpx2[ii,  jj] = um*( coeff*( - ucat_x[ii, jj+2] - 2*ucat_x[ii, jj+1] + 3*ucat_x[ii, jj] ) + ucat_x[ii, jj+1] ) + up*( coeff*( - ucat_x[ii, jj-1] - 2*ucat_x[ii, jj] + 3*ucat_x[ii, jj+1] ) + ucat_x[ii, jj] )
                  
                  fpy2[ii,  jj] = vm*( coeff*( - ucat_y[ii, jj+2] - 2*ucat_y[ii, jj+1] + 3*ucat_y[ii, jj] ) + ucat_y[ii, jj+1] ) + vp*( coeff*( - ucat_y[ii, jj-1] - 2*ucat_y[ii, jj] + 3*ucat_y[ii, jj+1] ) + ucat_y[ii, jj] )
               end
            end
         end

         for ii = 1:M
            fpx2[ii, N] = 0
            fpy2[ii, N] = 0
         end

         for ii = 2:M-1
            for jj = 2:N-1
               conv_flux_x[ii, jj] = ( fpx1[ii, jj] - fpx1[ii-1, jj] )/dx + ( fpx2[ii, jj] - fpx2[ii, jj-1] )/dy
               conv_flux_y[ii, jj] = ( fpy1[ii, jj] - fpy1[ii-1, jj] )/dx + ( fpy2[ii, jj] - fpy2[ii, jj-1] )/dy
            end
         end

         for ii = 1:M
            conv_flux_x[ii, 1] = 0.0
            conv_flux_y[ii, 1] = 0.0
            conv_flux_x[ii, N] = 0.0
            conv_flux_y[ii, N] = 0.0
         end

         for jj = 1:N
            conv_flux_x[1, jj] = 0.0
            conv_flux_y[1, jj] = 0.0
            conv_flux_x[M, jj] = 0.0
            conv_flux_y[M, jj] = 0.0
         end

         return conv_flux_x, conv_flux_y
      end
   end
end