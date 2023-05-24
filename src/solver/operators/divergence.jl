module Operators

export divergence

function divergence(ucurx, ucury, dx, dy, M, N)
   div = Matrix{Float64}(undef, M, N)
   ucur_x = ucurx
   ucur_y = ucury

   for ii = 2:M-1
      for jj = 2:N-1
         div[ii, jj] = ( ucur_x[ii, jj] - ucur_x[ii-1, jj] )/dx + ( ucur_y[ii, jj] - ucur_y[ii, jj-1] )/dy
      end
   end
  
   # Zero out boundaries
   for runi = 1:M
      div[runi, 1] = 0.0
      div[runi, N] = 0.0
   end

   for runj = 1:N
      div[1, runj] = 0.0
      div[M, runj] = 0.0
   end

   return div
end

end