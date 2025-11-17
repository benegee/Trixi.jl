
using Statistics
using BenchmarkTools
using OrdinaryDiffEq
using Trixi
using CUDA
CUDA.allowscalar(false)

function error_statistics(reference, actual)
      diff = reference - actual
      println("\tMin error:    ", minimum(diff))
      println("\tMax error:    ", maximum(diff))
      println("\tMean error:   ", Statistics.mean(diff))
      println("\tMedian error: ", Statistics.median(diff))
end

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations3D(1.4)

function initial_condition_taylor_green_vortex(x, t,
                                               equations::CompressibleEulerEquations3D)
    A  = 1.0 # magnitude of speed
    Ms = 0.1 # maximum Mach number

    rho = 1.0
    v1  =  A * sin(x[1]) * cos(x[2]) * cos(x[3])
    v2  = -A * cos(x[1]) * sin(x[2]) * cos(x[3])
    v3  = 0.0
    p   = (A / Ms)^2 * rho / equations.gamma # scaling to get Ms
    p   = p + 1.0/16.0 * A^2 * rho * (cos(2*x[1])*cos(2*x[3]) +
          2*cos(2*x[2]) + 2*cos(2*x[1]) + cos(2*x[2])*cos(2*x[3]))

    return prim2cons(SVector(rho, v1, v2, v3, p), equations)
end

initial_condition = initial_condition_taylor_green_vortex
surface_flux = flux_lax_friedrichs
volume_flux = flux_kennedy_gruber
solver = DGSEM(polydeg=3, surface_flux=surface_flux,
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (-1.0, -1.0, -1.0) .* pi
coordinates_max = ( 1.0,  1.0,  1.0) .* pi

initial_refinement_level = 1
trees_per_dimension = (4, 4, 4)

mesh = P4estMesh(trees_per_dimension, polydeg=3,
                 coordinates_min=coordinates_min, coordinates_max=coordinates_max,
                 periodicity=true, initial_refinement_level=initial_refinement_level)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan; adapt_to=CuArray)

summary_callback = SummaryCallback()

stepsize_callback = StepsizeCallback(cfl=0.1)

callbacks = CallbackSet(summary_callback,
                        stepsize_callback)


###############################################################################
# run the simulation

integrator = init(ode, CarpenterKennedy2N54(williamson_condition=false),
                  dt=1.0,
                  save_everystep=false, callback=callbacks)
solve!(integrator)

mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(ode.p)

u_ode = integrator.u
u = Trixi.wrap_array(u_ode, mesh, equations, solver, cache)
du_ref = similar(u)
du_ref .= 0
du_new = similar(u)
du_new .= 0

Trixi.calc_volume_integral!(du_ref, u, mesh, Trixi.False(), equations, solver.volume_integral, solver, cache)
Trixi.calc_volume_integral!(du_new, u, mesh, Trixi.False(), equations, solver.volume_integral, solver, cache)

if all(du_ref .≈ du_new)
      println("Sanity check passed.")
else
      println("[ERR] Sanity check FAILED.")
      error_statistics(du_ref, du_new)
end

# Testing
println()
du_exp = similar(u)
du_exp .= 0
Trixi.exp_index_calc_volume_integral!(du_exp, u, mesh, Trixi.False(), equations, solver.volume_integral, solver, cache)

if all(du_ref .≈ du_exp)
      println("The exp_index version is corrrect.")
else
      println("[ERR] There is a BUG in the exp_index version.")
      error_statistics(du_ref, du_exp)
end

du_exp = similar(u)
du_exp .= 0
Trixi.exp_ijk_calc_volume_integral!(du_exp, u, mesh, Trixi.False(), equations, solver.volume_integral, solver, cache)

if all(du_ref .≈ du_exp)
      println("The exp_ijk version is corrrect.")
else
      println("[ERR] There is a BUG in the exp_ijk version.")
      error_statistics(du_ref, du_exp)
end

# Tuning
println()
best_time = Inf
wgs = 0
best_wgs = wgs
index_x = 1
println("Tuning reference")
while index_x * 32 <= 1024
      global wgs = index_x * 32
      try
            res = @btimed begin
                  Trixi.calc_volume_integral!(du_ref, u, mesh, Trixi.False(), equations, solver.volume_integral, solver, cache, wgs)
                  CUDA.device_synchronize()
            end
            if res.time < best_time
                  global best_time = res.time
                  global best_wgs = wgs
            end
            global index_x += 1
      catch
            global index_x += 1
      end
end
println("\tBest time: ", best_time, " s -- workgroupsize: ", best_wgs)
println("Tuning exp_index")
best_time = Inf
wgs = 0
best_wgs = wgs
index_x = 1
while index_x * 32 <= 1024
      global wgs = index_x * 32
      try
            res = @btimed begin
                  Trixi.exp_index_calc_volume_integral!(du_exp, u, mesh, Trixi.False(), equations, solver.volume_integral, solver, cache, wgs)
                  CUDA.device_synchronize()
            end
            if res.time < best_time
                  global best_time = res.time
                  global best_wgs = wgs
            end
            global index_x += 1
      catch
            global index_x += 1
      end
end
println("\tBest time: ", best_time, " s -- workgroupsize: ", best_wgs)
println("Tuning exp_ijk")
best_time = Inf
wgs = (0, 0)
best_wgs = wgs
index_x = 1
index_y = 1
while index_x * 32 <= 1024
      while index_y <= 32
            if index_x * index_y > 1024
                  global index_x += 1
                  global index_y = 1
                  continue
            end
            global wgs = (index_x * 32, index_y)
            try
                  res = @btimed begin
                        Trixi.exp_ijk_calc_volume_integral!(du_exp, u, mesh, Trixi.False(), equations, solver.volume_integral, solver, cache, wgs)
                        CUDA.device_synchronize()
                  end
                  if res.time < best_time
                        global best_time = res.time
                        global best_wgs = wgs
                  end
                  global index_y += 1
            catch
                  global index_y += 1
            end
      end
      global index_x += 1
      global index_y = 1
end
println("\tBest time: ", best_time, " s -- workgroupsize: ", best_wgs)

finalize(mesh)
