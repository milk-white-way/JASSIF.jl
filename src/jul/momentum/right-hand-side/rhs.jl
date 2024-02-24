module RHS

include("./convective_flux.jl")
include("./viscous_flux.jl")
include("./pressure_grad.jl")

using LinearAlgebra
using .ConvectiveFlux
using .ViscousFlux
using .PressureGradient
 
export rhs_calc

	function rhs_newton(CSolution, convect_coeff)
	   form_bcs(CSolution)
      rhs_calc(CSolution)

      M = CSolution.m2
      N = CSolution.n2

      dt = CSolution.dt
      local_rhs_x = CSolution.rhs_x
      local_rhs_y = CSolution.rhs_y

      local_ucur_x = CSolution.ucur_x
      local_ucur_y = CSolution.ucur_y

      local_rhs_x = local_rhs_x - ( 1.5/dt )*( uim_x - local_ucur_x ) + ( 0.5/dt )*(  )
      local_rhs_y = local_rhs_y - ( 1.5/dt )*( uim_y - local_ucur_y ) + ( 0.5/dt )*(  )

      for ii = 1:M
	      local_rhs_y[ii,   1] = 0
	      local_rhs_y[ii,   N] = 0
	      local_rhs_y[ii, N-1] = 0
      end

      for jj = 1:N
	      local_rhs_x[  1, jj] = 0
	      local_rhs_x[  M, jj] = 0
	      local_rhs_x[M-1, jj] = 0
      end

      global CSolution.rhs_x = local_rhs_x
      global CSolution.rhs_y = local_rhs_y
   end

	function rhs_calc(CSolution, convec_coeff)
		#println("Calculating the Right Hand Side...")
		let M = deepcopy(CSolution.m2)
          N = deepcopy(CSolution.n2)
			 rhs_x = Matrix{Float64}(undef, M, N)
			 rhs_y = Matrix{Float64}(undef, M, N)

			conv_flux_x, conv_flux_y = rhs_convective_flux_calc(CSolution, convec_coeff)
			visc_flux_x, visc_flux_y = rhs_viscous_flux_calc(CSolution)
			pres_grad_x, pres_grad_y = rhs_pressure_gradient_calc(CSolution)

			total_flux_x = - conv_flux_x + visc_flux_x - pres_grad_x
			total_flux_y = - conv_flux_y + visc_flux_y - pres_grad_y

			for ii = 1:M-1
				for jj = 1:N
					rhs_x[ii, jj] = 0.5*( total_flux_x[ii, jj] + total_flux_x[ii+1, jj] )
				end
			end

			for ii = 1:M
				for jj = 1:N-1
					rhs_y[ii, jj] = 0.5*( total_flux_y[ii, jj] + total_flux_y[ii, jj+1] )
				end
			end

			for ii = 1:M
				rhs_x[ii,   1] = 0
				rhs_x[ii,   N] = 0

				rhs_y[ii,   1] = 0
				rhs_y[ii, N-1] = 0
				rhs_y[ii,   N] = 0
			end

			for jj = 1:N
				rhs_x[  M, jj] = 0
				rhs_x[M-1, jj] = 0
				rhs_x[  1, jj] = 0

				rhs_y[  1, jj] = 0
				rhs_y[  M, jj] = 0
			end

			return rhs_x, rhs_y
		end
		#println("Righ Hand Side Calculated!")
	end
end
