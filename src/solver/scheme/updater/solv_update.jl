module Updater

include("../../operators/lidx.jl")

using .LIDX

export solv_update

function solv_update(CSolution, phi)
   let M = deepcopy(CSolution.m2)
       N = deepcopy(CSolution.n2)
       ucur_x = deepcopy(CSolution.ucur_x)
       ucur_y = deepcopy(CSolution.ucur_y)
       pres = deepcopy(CSolution.pres)
       dx = deepcopy(CSolution.dx)
       dy = deepcopy(CSolution.dy)
       dt = deepcopy(CSolution.dt)
  
      a_Phi = Matrix{Float64}(undef, M, N)
      P_G_x = Matrix{Float64}(undef, M, N)
      P_G_y = Matrix{Float64}(undef, M, N)

      for ii = 1:M
         for jj = 1:N
            a_Phi[ii, jj] = 0.0
            P_G_x[ii, jj] = 0.0
            P_G_y[ii, jj] = 0.0
         end
      end

      for counter = 1:M*N
         ii, jj = lidx(counter, M, N)
         a_Phi[ii, jj] = phi[counter]
      end

      for ii = 2:M-2
         for jj = 2:N-1
            P_G_x[ii, jj] = ( a_Phi[ii+1, jj] - a_Phi[ii, jj] )/dx
         end
      end
  
      for ii = 2:M-1
         for jj = 2:N-2
            P_G_y[ii, jj] =  ( a_Phi[ii, jj+1] - a_Phi[ii, jj] )/dy
         end
      end
  
      ucur_x_new = deepcopy(ucur_x)
      ucur_y_new = deepcopy(ucur_y)

      for ii = 1:M-1
         for jj = 1:N
            ucur_x_new[ii, jj] = ucur_x[ii, jj] - P_G_x[ii, jj]*dt/1.5
         end
      end
  
      for ii = 1:M
         for jj = 1:N-1
            ucur_y_new[ii, jj] = ucur_y[ii, jj] - P_G_y[ii, jj]*dt/1.5
         end
      end
  
      pres_new = pres + a_Phi
  
      global CSolution.ucur_x = ucur_x_new
      global CSolution.ucur_y = ucur_y_new
      global CSolution.pres = pres_new
   end
end

end