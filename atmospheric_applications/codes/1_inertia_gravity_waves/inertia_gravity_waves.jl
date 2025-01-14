using OrdinaryDiffEq
using Trixi
using Plots

# Definition of struct field to collect all the physical parameters and hydrostatic background state

struct InertiaSetup
    # Physical constants
    g::Float64       # gravity of earth
    c_p::Float64     # heat capacity for constant pressure (dry air)
    c_v::Float64     # heat capacity for constant volume (dry air)
    gamma::Float64   # heat capacity ratio (dry air)
    p_0::Float64     # atmospheric pressure
    theta_0::Float64  # constant background potential temperature steady state     
    function InertiaSetup(; g = 9.81, c_p = 1004.0, c_v = 717.0, gamma = c_p / c_v, p_0 = 100_000.0, theta_0 = 300.0)
        new(g, c_p, c_v, gamma, p_0, theta_0)
    end
end

# Initial condition for inertia gravity waves: constant potential temperature background state plus a perturbation
function (setup::InertiaSetup)(x, t, equations::CompressibleEulerEquations2D)
    
    @unpack g, c_p, c_v, p_0, theta_0 = setup
    
    # center of perturbation
    x_c = 100_000.0
    a = 5_000
    H = 10_000
    
    # perturbation in potential temperature
    potential_temperature_ref = theta_0 * exp(0.01^2 / g * x[2])
    potential_temperature_perturbation = 0.01 * sinpi(x[2] / H) / (1 + (x[1] - x_c)^2 / a^2)
    potential_temperature = potential_temperature_ref + potential_temperature_perturbation

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 + g^2 / (c_p * 300.0 * 0.01^2) * (exp(-0.01^2 / g * x[2]) - 1)

    # pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)
    T = potential_temperature * exner

    # density
    rho = p / (R * T)
    v1 = 20.0
    v2 = 0.0

    return prim2cons(SVector(rho, v1, v2, p), equations)
end

# Source terms for the CompressibleEulerEquations2D
@inline function (setup::InertiaSetup)(u, x, t, equations::CompressibleEulerEquations2D)
    @unpack g = setup
    rho, _, rho_v2, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -g * rho, -g * rho_v2)
end

###############################################################################
# Semidiscretization of the compressible Euler equations
inertia_gravity_wave_setup = InertiaSetup()

equations = CompressibleEulerEquations2D(inertia_gravity_wave_setup.gamma)

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = FluxLMARS(340.0)

volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)

cells_per_dimension = (32, 32)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))

semi = SemidiscretizationHyperbolic(mesh, equations, inertia_gravity_wave_setup, solver,
                                    source_terms = inertia_gravity_wave_setup,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 3000.0)  # 1000 seconds final time

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "1_inertia_gravity_waves/out",
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation
sol = solve(ode, SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);

summary_callback()

# Visulization: 

# Compute the potential temperature perturbation
vertical_velocity = let u = Trixi.wrap_array(sol.u[end], semi)
    @unpack g, c_p, c_v, p_0, theta_0 = inertia_gravity_wave_setup
    rho = u[1, :, : ,:]
    rho_v2 = u[3, : ,: ,:]
    rho_v2 ./ rho
end

# Plotting the vertical velocity component
plot(ScalarPlotData2D(vertical_velocity, semi), title = "Vertical velocity component [m/s]", aspect_ratio = 5)



