module Continuity

include("../operators/glidx.jl")
include("../operators/lidx.jl")
include("../operators/divergence.jl")

using TickTock
using LinearAlgebra
using SparseArrays
using Printf
using Debugger
using .GLIDX
using .LIDX
using .Operators

function poisson_solv(CSolution)
   let M = deepcopy(CSolution.m2)
       N = deepcopy(CSolution.n2)
       ucur_x = deepcopy(CSolution.ucur_x)
       ucur_y = deepcopy(CSolution.ucur_y)
       dx = deepcopy(CSolution.dx)
       dy = deepcopy(CSolution.dy)
       dt = deepcopy(CSolution.dt)

      A = poisson_lhs_neumann(M, N, dx, dy)
      #println(spy(A))
      b = poisson_rhs_newmann(ucur_x, ucur_y, dx, dy, dt, M, N)
      phi = A\b'
      return phi
  end
end

function poisson_lhs_neumann(M, N, dx, dy)
   let A = spzeros(M*N, M*N)
       sp_ii = Float64[]
       sp_jj = Float64[]
       sp_vals = Float64[]

      for counter in 1:M*N
         #tick()
         ii, jj = lidx(counter, M, N)

         indexN = glidx(ii  , jj-1, M, N)
         indexE = glidx(ii-1, jj  , M, N)
         indexP = glidx(ii  , jj  , M, N)
         indexW = glidx(ii+1, jj  , M, N)
         indexS = glidx(ii  , jj+1, M, N)
      
         if ( ( ii>=2 ) && ( ii<=M-1 ) && ( jj>=2 ) && ( jj<=N-1 ) )
            #ii = []
            #jj = []
            #vals = []

            push!( sp_ii, indexP );push!( sp_ii, indexP );push!( sp_ii, indexP );push!( sp_ii, indexP );push!( sp_ii, indexP )
            push!( sp_jj, indexN );push!( sp_jj, indexE );push!( sp_jj, indexP );push!( sp_jj, indexW );push!( sp_jj, indexS )
            push!( sp_vals, 1/( dy^2 ) );push!( sp_vals, 1/( dx^2 ) );push!( sp_vals, -2/( dx^2 ) - 2/( dy^2 ) );push!( sp_vals, 1/( dx^2 ) );push!( sp_vals, 1/( dy^2 ) )
            
            if ( ( ii==2 ) && ( jj==2 ) )
               sp_vals[3] = 0
            end

         else
            if ( ( ii==1 ) && ( jj!=1 ) && ( jj!=N ) )
               #ii = []
               #jj = []
               #vals = []
        
               push!( sp_ii, indexP );push!( sp_ii, indexP )
               push!( sp_jj, indexP );push!( sp_jj, indexW )
               push!( sp_vals, -1 ); push!( sp_vals, 1 )
            end
      
            if ( ( ii==M ) && ( jj!=1 ) && ( jj!=N ) )
               #ii = []
               #jj = []
               #vals = []
        
               push!( sp_ii, indexP );push!( sp_ii, indexP )
               push!( sp_jj, indexP );push!( sp_jj, indexE )
               push!( sp_vals, -1 ); push!( sp_vals, 1 )
            end
      
            if ( ( jj==1 ) && ( ii!=1 ) && ( ii!=M ) )
               #ii = []
               #jj = []
               #vals = []
        
               push!( sp_ii, indexP );push!( sp_ii, indexP )
               push!( sp_jj, indexP );push!( sp_jj, indexS )
               push!( sp_vals, -1 ); push!( sp_vals, 1 )
            end
      
            if ( ( jj==N ) && ( ii!=1 ) && ( ii!=M ) )
               #ii = []
               #jj = []
               #vals = []
        
               push!( sp_ii, indexP );push!( sp_ii, indexP )
               push!( sp_jj, indexP );push!( sp_jj, indexN )
               push!( sp_vals, -1 ); push!( sp_vals, 1 )
            end
      
            if ( ( ii==1 && jj==1 ) || ( ii==1 && jj==N ) || ( ii==M && jj==1 ) || ( ii==M && jj==N ) )
               #ii = []
               #jj = []
               #vals = []
        
               push!( sp_ii, indexP )
               push!( sp_jj, indexP )
               push!( sp_vals, 1 )
            end
         end
         #tock()
         #=
         if ( counter==1 )
            sp_ii = ii
            sp_jj = jj
            sp_vals = vals
         else
            sp_ii = vcat(sp_ii, ii)
            sp_jj = vcat(sp_jj, jj)
            sp_vals = vcat(sp_vals, vals)
         end
         =#
      end
      A = sparse(sp_ii, sp_jj, sp_vals)
      return A
   end
end

function poisson_rhs_newmann(ucur_x, ucur_y, dx, dy, dt, M, N)
   P_RHS = zeros(1, M*N)

   P_Div = divergence(ucur_x, ucur_y, dx, dy, M, N)

   for ii = 1:M
      for jj = 1:N
         index = glidx(ii, jj, M, N)
         if ( ( ii>=2 ) && ( ii<=M-1 ) && ( jj>=2 ) && ( jj<=N-1 ) )
            P_RHS[index] = P_Div[ii, jj]
         else
            # Apply Neumann boundary condition for Poisson equation
            P_RHS[index] = 0.0
         end
      end
   end
  
      # Check if solvability condition is met !
      Summation = 0
      for ii = 1:M
         for jj = 1:N
            index = glidx(ii, jj, M, N)
            if ( ( ii>=2 ) && ( ii<=M-1 ) && ( jj>=2 ) && ( jj<=N-1 ) )
               # Not the fix pressure location
               Summation = Summation + P_RHS[index]
            end
         end
      end
      Summation = ( 1.5*Summation )/dt
      if ( Summation > 1e-8 )
         @printf("Solvability condition not met ! Summation = %8.5f \n", Summation)
      end
  
      # Scale back to time discretization scheme
      # Adam-Bashforth scheme
      P_RHS = ( 1.5*P_RHS )/dt
      return P_RHS
end

end