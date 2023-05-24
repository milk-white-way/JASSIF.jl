module NewtonKrylov

using LinearAlgebra
using SparseArrays

export jfnk

function jfnk(CSolution, maxiter, toleran, convect_coeff)
  
  M = CSolution.m2
  N = CSolution.n2

  A = sparse(1.0I, 2*M*N, 2*M*N)

  zero_guess = zeros(2*M*N)

  # Input zero guess into RHS function
  # Calculate the value of RHS
  # Calulate norm
end

end