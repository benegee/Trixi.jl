using OrdinaryDiffEq
using Trixi
using CUDA
CUDA.allowscalar(false)

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
      println("The experimental version is corrrect.")
else
      println("There is a bug in the experimental version.")
end

finalize(mesh)