module Initialization

export init

function init(CSolution, SYS2D_LENGTH_X, SYS2D_LENGTH_Y,
               SYS2D_GHOST_X, SYS2D_GHOST_Y,
               SYS2D_STEP_X, SYS2D_STEP_Y,
               STEP_TIME, ENDO_TIME,
               REN_NUM, DYN_VIS)

   Ucat_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
   Ucat_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
   Ucur_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
   Ucur_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
   Ubcs_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
   Ubcs_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
   Uint_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
   Uint_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
     Pres = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
     dU_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
     dU_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)

   for ii = 1:SYS2D_GHOST_X
      for jj = 1:SYS2D_GHOST_Y
         Ucat_x[ii, jj] = 0
         Ucat_y[ii, jj] = 0
         Ucur_x[ii, jj] = 0
         Ucur_y[ii, jj] = 0
         Ubcs_x[ii, jj] = 0
         Ubcs_y[ii, jj] = 0
         Uint_x[ii, jj] = 0
         Uint_y[ii, jj] = 0
           Pres[ii, jj] = 0
           dU_x[ii, jj] = 0
           dU_y[ii, jj] = 0
      end
   end

   global CSolution.length_x = SYS2D_LENGTH_X
   global CSolution.length_y = SYS2D_LENGTH_Y
   global CSolution.m2 = SYS2D_GHOST_X
   global CSolution.n2 = SYS2D_GHOST_Y
   global CSolution.dx = SYS2D_STEP_X
   global CSolution.dy = SYS2D_STEP_Y
   global CSolution.dt = STEP_TIME
   global CSolution.et = ENDO_TIME
   global CSolution.ren = REN_NUM
   global CSolution.vis = DYN_VIS
   global CSolution.ucat_x = Ucat_x
   global CSolution.ucat_y = Ucat_y
   global CSolution.ucur_x = Ucur_x
   global CSolution.ucur_y = Ucur_y
   global CSolution.ubcs_x = Ubcs_x
   global CSolution.ubcs_y = Ubcs_y
   global CSolution.uint_x = Uint_x
   global CSolution.uint_y = Uint_y
   global CSolution.pres = Pres
   global CSolution.du_x = dU_x
   global CSolution.du_y = dU_y
end

end
