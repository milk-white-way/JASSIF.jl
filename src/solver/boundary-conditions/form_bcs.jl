module FormBCS

   include("../operators/cur2cat.jl")
   using .Converts
   using Printf

   export form_bcs

   function form_bcs(CSolution)

      let M = deepcopy(CSolution.m2)
          N = deepcopy(CSolution.n2)
          dx = deepcopy(CSolution.dx)
          dy = deepcopy(CSolution.dy)
          ucur_x = deepcopy(CSolution.ucur_x)
          ucur_y = deepcopy(CSolution.ucur_y)

          ubcs_x_1 = zeros(Float64, 1, M)
          ubcs_x_2 = zeros(Float64, 1, M)
          ubcs_y_1 = zeros(Float64, 1, M)
          ubcs_y_2 = zeros(Float64, 1, M)

          ubcs_x_3 = zeros(Float64, 1, N)
          ubcs_x_4 = zeros(Float64, 1, M)
          ubcs_y_3 = zeros(Float64, 1, M)
          ubcs_y_4 = zeros(Float64, 1, M)
   
         for jj = 1:N
            ucur_x[  1, jj] = 0.0
            ucur_x[M-1, jj] = 0.0
         end
       
         for ii = 1:M
            ucur_y[ii,   1] = 0.0
            ucur_y[ii, N-1] = 0.0
         end

         for jj = 1:N
            ii = 1
            x = ( ii-2+0.5 )*dx
            y = ( jj-2 )*dy
            ubcs_x_3[jj] = 0.0
            ubcs_y_3[jj] = 0.0

            ii = M-1
            x = (jj-2+0.5)*dx
            y = (jj-2 )*dy
            ubcs_x_4[jj] = 0.0
            ubcs_y_4[jj] = 0.0
         end

         for ii = 1:M
            jj = 1
            x = ( ii-2 )*dx
            y = ( jj-2+0.5 )*dy
            ubcs_x_1[ii] = 0.0
            ubcs_y_1[ii] = 0.0

            jj = N-1
            x = ( ii-2 )*dx
            y = ( jj-2+0.5 )*dy
            ubcs_x_2[ii] = 1.0
            ubcs_y_2[ii] = 0.0

         end

         #println(ubcs_x_2)
         global CSolution.ucur_x = ucur_x
         global CSolution.ucur_y = ucur_y

         cur2cat(CSolution)

         ucat_x = deepcopy(CSolution.ucat_x)
         ucat_y = deepcopy(CSolution.ucat_y)
  
         for ii = 1:M
            ucat_x[ii, 1] = 2*ubcs_x_1[ii] - ucat_x[ii, 2]
            ucat_y[ii, 1] = 2*ubcs_y_1[ii] - ucat_y[ii, 2]

            ucat_x[ii, N] = 2*ubcs_x_2[ii] - ucat_x[ii, N-1]
            ucat_y[ii, N] = 2*ubcs_y_2[ii] - ucat_y[ii, N-1]
         end

         for jj = 1:N
            ucat_x[1, jj] = 2*ubcs_x_3[jj] - ucat_x[2, jj]
            ucat_y[1, jj] = 2*ubcs_y_3[jj] - ucat_y[2, jj]

            ucat_x[M, jj] = 2*ubcs_x_4[jj] - ucat_x[M-1, jj]
            ucat_y[M, jj] = 2*ubcs_y_4[jj] - ucat_y[M-1, jj]
         end

         ucat_x[1, 1] = 0.0
         ucat_y[1, 1] = 0.0

         ucat_x[1, N] = 0.0
         ucat_y[1, N] = 0.0

         ucat_x[M, 1] = 0.0
         ucat_y[M, 1] = 0.0

         ucat_x[M, N] = 0.0
         ucat_y[M, N] = 0.0

         global CSolution.ucat_x = ucat_x
         global CSolution.ucat_y = ucat_y
      end
   end
end
