# A turbulent helically driven dynamo model in resistive compressible MHD.
# We start with a weak small-scale magnetic field provided as Gaussian noise.
# The energy is supplied via a forcing function f_force(x, t), which is delta-correlated in time.
# For a comparison look into the Pencil Code example https://github.com/pencil-code/pencil-code/tree/master/samples/helical-MHDturb.
# For more information on the initial condition and forcing term see
#  S. Candelaresi and A. Brandenburg (2013)
#  Kinetic helicity needed to drive large-scale dynamos
#  [DOI: 10.1103/PhysRevE.87.043104](https://doi.org/10.1103/PhysRevE.87.043104)
#  section II A

using OrdinaryDiffEq
using Trixi
using Random
using LinearAlgebra
using StableRNGs

###############################################################################
# initial condition and forcing function

"""
    initial_condition_gaussian_noise_mhd(x, t, equations)

Define the initial condition with a weak small-scale Gaussian noise magnetic field.
"""
function initial_condition_gaussian_noise_mhd(x, t, equations)
    # Random weak small-scale initial magnetic field.
    amplitude = 1e-8

    seed = reinterpret(UInt, t + x[1] + x[2] + x[3])
    rng = StableRNG(seed)

    rho = 1.0
    rho_v1 = 0.0
    rho_v2 = 0.0
    rho_v3 = 0.0
    rho_e = 10.0
    B1 = randn(rng) * amplitude
    B2 = randn(rng) * amplitude
    B3 = randn(rng) * amplitude
    psi = 0.0

    return SVector(rho, rho_v1, rho_v2, rho_v3, rho_e, B1, B2, B3, psi)
end

"""
    source_terms_helical_forcing(u, x, t, equations::ViscoResistiveMhd3D)

Forcing term that adds a helical small-scale driver to the system that is
delta-correlated in time.
For more information on the initial condition and forcing term see 
- S. Candelaresi and A. Brandenburg (2013) 
  Kinetic helicity needed to drive large-scale dynamos
  [DOI: 10.1103/PhysRevE.87.043104](https://doi.org/10.1103/PhysRevE.87.043104)
"""
function source_terms_helical_forcing(u, x, t, equations::IdealGlmMhdEquations3D)
    # forcing amplitude
    f_0 = 0.02
    # helicality of the forcing (-1, 1)
    sigma = 1.0

    # Extract some parameters for the computation.
    c_s = 0.1
    delta_t = 0.01

    # To make sure that the random numbers depend only on time we need to set the seeds.
    seed = reinterpret(UInt, t)
    rng = StableRNG(seed)

    # Random phase -pi < phi <= pi
    phi = (randn(rng) * 2 - 1) * pi

    # Random vector k, also delta-correlated in time.
    k = SVector((randn(rng) * 2 - 1) * pi, (randn(rng) * 2 - 1) * pi, (randn(rng) * 2 - 1) * pi)
    k = k / norm(k)
    k = 6 * k
    k_norm = k / norm(k)

    # Random unit vector, not aligned with k.
    ee = SVector(randn(rng), randn(rng), randn(rng)) .* 2 .- 1
    # Normalize the vector e.
    ee = ee / norm(ee)

    # Compute the f-vectors.
    f_nohel = cross(k, ee) / sqrt(dot(k, k) - dot(k, ee)^2)
    M = [0 k_norm[3] -k_norm[2]; -k_norm[3] 0 k_norm[1]; k_norm[2] -k_norm[1] 0]
    R = (I - im * sigma * M) / sqrt(1 + sigma^2)
    f_k = R * f_nohel

    # Normalization factor to make sure that the time integration of f is 0.
    N = f_0 * c_s * sqrt(norm(k) * c_s / delta_t)

    forcing = real(N * f_k * exp(im * dot(k, x) + im * phi)) / u[1]

    return SVector(0, forcing[1], forcing[2], forcing[3], 0, 0, 0, 0, 0)
end


###############################################################################
# semidiscretization of the visco-resistive compressible MHD equations

prandtl_number() = 0.72
mu() = 5e-3
eta() = 5e-3

# Adiabatic monatomic gas in 3d 
gamma = 1.0 + 2.0 / 3.0

equations = IdealGlmMhdEquations3D(gamma)
equations_parabolic = ViscoResistiveMhd3D(equations, mu = mu(),
					  Prandtl = prandtl_number(), eta = eta(),
                                          gradient_variables = GradientVariablesPrimitive())

volume_flux = (flux_hindenlang_gassner, flux_nonconservative_powell)
solver = DGSEM(polydeg = 3,
               surface_flux = (flux_lax_friedrichs, flux_nonconservative_powell),
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (-pi / 4, -pi / 4, 0.0) # minimum coordinates (min(x), min(y), min(z))
coordinates_max = (pi / 4, pi / 4, pi / 2) # maximum coordinates (max(x), max(y), max(z))

# Create a uniformly refined mesh with periodic boundaries
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                periodicity = (true, true, true),
                n_cells_max = 100_000) # set maximum capacity of tree data structure

initial_condition = initial_condition_gaussian_noise_mhd
source_terms = source_terms_helical_forcing

semi = SemidiscretizationHyperbolicParabolic(mesh, (equations, equations_parabolic),
                                             initial_condition, solver,
                                             source_terms = source_terms)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span `tspan`
tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
alive_callback = AliveCallback(alive_interval = 10)
analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)
save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

save_restart = SaveRestartCallback(interval = 100,
                                   save_final_restart = true)

cfl = 1.5
stepsize_callback = StepsizeCallback(cfl = cfl)

glm_speed_callback = GlmSpeedCallback(glm_scale = 0.5, cfl = cfl)

callbacks = CallbackSet(summary_callback,
                        alive_callback,
                        analysis_callback,
                        save_solution,
                        stepsize_callback,
                        glm_speed_callback,
                        save_restart)

###############################################################################
# run the simulation

time_int_tol = 1e-5
sol = solve(ode, RDPK3SpFSAL49(), dt = 1e-5, save_everystep = false, callback = callbacks)
summary_callback() # print the timer summary