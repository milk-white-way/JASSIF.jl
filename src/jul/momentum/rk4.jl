module Momentum

   include("../boundary-conditions/form_bcs.jl")
   include("../momentum/right-hand-side/rhs.jl")

   using Printf
   using LinearAlgebra
   using Debugger # Not being used
   using .FormBCS
   using .RHS

   export runge_kutta

   function runge_kutta(CSolution, maxiter, tol, convec_coeff)
      let coefficient = [1/4 1/3 1/2 1]
          ucur_x = deepcopy(CSolution.ucur_x)
          ucur_y = deepcopy(CSolution.ucur_y)
          du_x = deepcopy(CSolution.du_x)
          du_y = deepcopy(CSolution.du_y)
          u_p_x = deepcopy(ucur_x)
          u_p_y = deepcopy(ucur_y)

          iter_stp = 0
          iter_err = 1

         while ( iter_stp < maxiter ) && ( iter_err > tol )
            u_im_x = deepcopy(u_p_x)
            u_im_y = deepcopy(u_p_y)
            dt = deepcopy(CSolution.dt)

            for stage in 1:4
               rhs_x, rhs_y = rhs_calc(CSolution, convec_coeff)

               rhs_x = rhs_x - ( 1.5/dt )*( u_im_x - ucur_x ) + ( 0.5/dt )*( du_x )
               rhs_y = rhs_y - ( 1.5/dt )*( u_im_y - ucur_y ) + ( 0.5/dt )*( du_y )

               # Update Immidiate Velocity components
               u_im_x = u_p_x + coefficient[stage]*dt*0.4*rhs_x
               u_im_y = u_p_y + coefficient[stage]*dt*0.4*rhs_y

               global CSolution.ucur_x = u_im_x
               global CSolution.ucur_y = u_im_y

               form_bcs(CSolution)
            end
            iter_err = norm(u_p_x - u_im_x, Inf)
            iter_stp = iter_stp + 1
    
            u_p_x = u_im_x
            u_p_y = u_im_y
    
            @printf("@subitr = %i, Momentum convergence = %e \n", iter_stp, iter_err)
         end
      end
   end
end