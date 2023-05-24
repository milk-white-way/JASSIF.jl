## JULIAN CODE BY TAM THIEN NGUYEN
# ++ email: tam.thien.nguyen@ndus.edu
## The goal is to numerically solve the Navier-Stokes equation
# Assumption:
# ++ Fluid Mechanics:
# ====== Incompressible
# ====== <Item 2>
# ++ Fluid Dynamics:
# ====== 2 Dimensional
# ====== <Item 1>
# ++ Input
# ====== <Add input here>
# ++ Output
# ====== u_cart_x: x-component of the Cartesian velocity
# ====== u_cart_y: y-component of the Cartesian velocity
# ====== pressure: scalar pressure field
# ====== time: time
# =========================================================
## GETTING LIBRARIES
#using Pkg
#Pkg.add("LinearAlgebra")
#Pkg.add("Debugger")
# Loading pre-built modules from libraries
using LinearAlgebra
using SparseArrays
using Printf
using TickTock
using DelimitedFiles
#using Debugger
#using Plots
# =========================================================
# Loading local modules
# 'Operators' component
include("./operators/divergence.jl")
include("./operators/init.jl")
include("./operators/cur2cat.jl")
# 'Boundary Conditions' component
include("./boundary-conditions/form_bcs.jl")
# 'Scheme' component
include("./scheme/updater/solv_update.jl")
# 'Momentum' component
#include("./momentum/jfnk.jl")
include("./momentum/rk4.jl")
include("./momentum/right-hand-side/rhs.jl")
# 'Continuity' component
include("./continuity/poisson.jl")
# ==========================================================

## SCHEME CONSTANTS
const NK_MAX_ITERATION = 100
const NK_TOLERANCE = 1.0e-5
const RK_MAX_ITERATION = 20
const RK_TOLERANCE = 1.0e-7

## PROBLEM CONSTANTS
const SYS2D_LENGTH_X = 1.0
const SYS2D_LENGTH_Y = 1.0
const SYS2D_DIM_X = 61 
const SYS2D_DIM_Y = 61
const SYS2D_GHOST_X = SYS2D_DIM_X + 2
const SYS2D_GHOST_Y = SYS2D_DIM_Y + 2
const SYS2D_STEP_X = SYS2D_LENGTH_X / (SYS2D_DIM_X - 1)
const SYS2D_STEP_Y = SYS2D_LENGTH_Y / (SYS2D_DIM_Y - 1)
const STEP_TIME = 0.008
const ENDO_TIME = 4000 # seconds
const REN_NUM = 100
const DYN_VIS = 1
const CONVECTIVE_COEFFICIENT = 0.128 # 1/8
const OUTPUT_FREQUENCY = 2000
# Creating solution structure containing:
#   Constants: lengths, spatial dimensions, spatial steps, time step, end time, reynolds number, viscosity
#   Variables: cartesian velocity field, curvilinear velocity field, boundary conditions, initial velocity field (in cartesian), pressure field
#mutable struct Constant2D
mutable struct Solution2D
   length_x::Float64
   length_y::Float64
   m2::UInt64
   n2::UInt64
   dx::Float64
   dy::Float64
   dt::Float64
   et::Float64
   ren::UInt64
   vis::Float64

   ucat_x::Matrix{Float64}
   ucat_y::Matrix{Float64}
   ucur_x::Matrix{Float64}
   ucur_y::Matrix{Float64}
   ubcs_x::Matrix{Float64}
   ubcs_y::Matrix{Float64}
   uint_x::Matrix{Float64}
   uint_y::Matrix{Float64}
     pres::Matrix{Float64}
     du_x::Matrix{Float64}
     du_y::Matrix{Float64}
end
#=
   conv_flux_x::Matrix{Float64}
   conv_flux_y::Matrix{Float64}
   visc_flux_x::Matrix{Float64}
   visc_flux_y::Matrix{Float64}
   pres_grad_x::Matrix{Float64}
   pres_grad_y::Matrix{Float64}
   rhs_x::Matrix{Float64}
   rhs_y::Matrix{Float64}
end
=#
tick()
CSolution = Solution2D( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y), 
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                        Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y) )

Initialization.init(CSolution,
                        SYS2D_LENGTH_X, SYS2D_LENGTH_Y,
                        SYS2D_GHOST_X, SYS2D_GHOST_Y,
                        SYS2D_STEP_X, SYS2D_STEP_Y,
                        STEP_TIME, ENDO_TIME,
                        REN_NUM, DYN_VIS)

global current_time::Float64  = 0.0
global export_alarm::Float64  = 0.0
global timestep::UInt64       = 1
global expostep::UInt64       = 1

while ( current_time < CSolution.et )
   @printf("Time step: %d \n", timestep)
   global current_time = timestep*CSolution.dt
   global export_alarm = expostep*OUTPUT_FREQUENCY
   let du_x = deepcopy(CSolution.du_x)
       du_y = deepcopy(CSolution.du_y)
       u_pre_x = deepcopy(CSolution.ucur_x)
       u_pre_y = deepcopy(CSolution.ucur_y)

      #jfnk(CSolution, MAX_ITERATION, TOLERANCE, CONVECTIVE_COEFFICIENT)
      Momentum.runge_kutta(CSolution, RK_MAX_ITERATION, RK_TOLERANCE, CONVECTIVE_COEFFICIENT)

      phi = Continuity.poisson_solv(CSolution)

      Updater.solv_update(CSolution, phi)

      u_new_x = deepcopy(CSolution.ucur_x)
      u_new_y = deepcopy(CSolution.ucur_y)

      global CSolution.du_x = u_new_x - u_pre_x
      global CSolution.du_y = u_new_y - u_pre_y
   end

   FormBCS.form_bcs(CSolution)

   #Div = Operators.divergence(CSolution)
   #MaxDiv = norm(Div, Inf)
   if timestep == export_alarm
      ufield_filename = string("ufield_at_", @sprintf("%d", export_alarm), ".csv")
      vfield_filename = string("vfield_at_", @sprintf("%d", export_alarm), ".csv")
      pfield_filename = string("pfield_at_", @sprintf("%d", export_alarm), ".csv")

      writedlm(ufield_filename, CSolution.ucat_x, ", ")
      writedlm(vfield_filename, CSolution.ucat_y, ", ")
      writedlm(pfield_filename, CSolution.pres, ", ")
      
      global expostep += 1
   end

   global timestep += 1
end
tock()
