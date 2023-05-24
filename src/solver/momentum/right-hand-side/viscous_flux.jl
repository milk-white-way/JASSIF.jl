module ViscousFlux
   export rhs_viscous_flux_calc

   function rhs_viscous_flux_calc(CSolution)
      let M = deepcopy(CSolution.m2)
          N = deepcopy(CSolution.n2)
          dx = deepcopy(CSolution.dx)
          dy = deepcopy(CSolution.dy)
          ucat_x = deepcopy(CSolution.ucat_x)
          ucat_y = deepcopy(CSolution.ucat_y)
          ucur_x = deepcopy(CSolution.ucur_x)
          ucur_y = deepcopy(CSolution.ucur_y)
          ren = deepcopy(CSolution.ren)

          visc_flux_x = Matrix{Float64}(undef, M, N)
          visc_flux_y = Matrix{Float64}(undef, M, N)

         for ii = 2:M-1
            for jj = 2:N-1
               ue = ucat_x[ii-1, jj]
               up = ucat_x[  ii, jj]
               uw = ucat_x[ii+1, jj]

               un = ucat_x[ii, jj-1]
               us = ucat_x[ii, jj+1]

               ve = ucat_y[ii-1, jj]
               vp = ucat_y[  ii, jj]
               vw = ucat_y[ii+1, jj]

               vn = ucat_y[ii, jj-1]
               vs = ucat_y[ii, jj+1]

               visc_flux_x[ii, jj] = ( uw - 2*up + ue )./( dx*dx ) + ( us - 2*up + un )./( dy*dy )
               visc_flux_y[ii, jj] = ( vw - 2*vp + ve )./( dx*dx ) + ( vs - 2*vp + vn )./( dy*dy )
            end
         end

         for ii = 1:M
            visc_flux_x[ii, 1] = 0
            visc_flux_y[ii, 1] = 0
            visc_flux_x[ii, N] = 0
            visc_flux_y[ii, N] = 0
         end
  
         for jj = 1:N
            visc_flux_x[1, jj] = 0
            visc_flux_y[1, jj] = 0
            visc_flux_x[1, jj] = 0
            visc_flux_y[1, jj] = 0
         end

         visc_flux_x = visc_flux_x/ren
         visc_flux_y = visc_flux_y/ren

         return visc_flux_x, visc_flux_y
      end
   end
end