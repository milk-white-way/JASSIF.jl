module Converts

export cur2cat

   function cur2cat(CSolution)
      let M = deepcopy(CSolution.m2)
          N = deepcopy(CSolution.n2)
          ucur_x = deepcopy(CSolution.ucur_x)
          ucur_y = deepcopy(CSolution.ucur_y)
          ucat_x = deepcopy(CSolution.ucat_x)
          ucat_y = deepcopy(CSolution.ucat_y)
         for ii = 2:M-1
            for jj = 2:N-1
               ucat_x[ii, jj] = 0.5*( ucur_x[ii, jj] + ucur_x[ii - 1, jj] )
               ucat_y[ii, jj] = 0.5*( ucur_y[ii, jj] + ucur_y[ii, jj - 1] )
            end
         end

         for ii = 1:M
            ucat_x[ii, 1] = 0.0
            ucat_x[ii, N] = 0.0

            ucat_y[ii, 1] = 0.0
            ucat_y[ii, N] = 0.0
         end

         for jj = 1:N
            ucat_x[1, jj] = 0.0
            ucat_x[M, jj] = 0.0

            ucat_y[1, jj] = 0.0
            ucat_y[M, jj] = 0.0
         end

         global CSolution.ucat_x = ucat_x
         global CSolution.ucat_y = ucat_y
      end
   end
end
