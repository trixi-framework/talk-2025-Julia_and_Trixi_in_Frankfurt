using OrdinaryDiffEq
using Trixi
using Plots

# Definition of struct field to collect all the physical parameters and hydrostatic background state

struct AtmoSetup
    # Physical constants
    g::Float64       # gravity of earth
    c_p::Float64     # heat capacity for constant pressure (dry air)
    c_v::Float64     # heat capacity for constant volume (dry air)
    gamma::Float64   # heat capacity ratio (dry air)
    p_0::Float64     # atmospheric pressure
    theta_0::Float64  # constant background potential temperature steady state     
    function AtmoSetup(; g = 9.81, c_p = 1004.0, c_v = 717.0, gamma = c_p / c_v, p_0 = 100_000.0, theta_0 = 300.0)
        new(g, c_p, c_v, gamma, p_0, theta_0)
    end
end

# Initial condition for density current: constant potential temperature background state plus a perturbation
function (setup::AtmoSetup)(x, t,
	equations::CompressibleEulerEquations2D)

    @unpack g, c_p, c_v, p_0, theta_0 = setup

    # center of perturbation
    center_x = 0.0
    center_z = 3000.0
    # radius of perturbation
    radius = 1.0
    x_r = 4000.0
    z_r = 2000.0

    # distance of current x to center of perturbation
    r = sqrt((x[1] - center_x)^2/x_r^2 + (x[2] - center_z)^2/z_r^2)

    # perturbation in potential temperature
    potential_temperature_ref = theta_0
    potential_temperature_perturbation = 0.0
    if r <= radius
        potential_temperature_perturbation = -15/2 * (1+cospi(r))
    end
    potential_temperature = potential_temperature_ref + potential_temperature_perturbation

    # Exner pressure, solves hydrostatic equation for z
    exner = 1 - g / (c_p * potential_temperature_ref) * x[2]

    # pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)
    T = potential_temperature * exner

    # density
    rho = p / (R * T)
    v1 = 0.0
    v2 = 0.0

	return prim2cons(SVector(rho, v1, v2, p), equations)
end

# Source terms for the CompressibleEulerEquations2D
@inline function (setup::AtmoSetup)(u, x, t, equations::CompressibleEulerEquations2D)
    @unpack g = setup
    rho, _, rho_v2, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -g * rho, -g * rho_v2)
end

@inline function cons2perturb(u, equations::CompressibleEulerEquations2D)
    @unpack c_v, c_p, p_0, theta_0 = inertia_gravity_wave_setup
    R = c_p - c_v
    rho, rho_v1, rho_v2, rho_e = u
    v1 = rho_v1 / rho
	v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
	rho_theta =  (p / p_0)^(c_v / c_p) * p_0 / R
    theta = rho_theta/rho - theta_0

    return  rho, v1, v2, theta

end

@inline Trixi.varnames(::typeof(cons2perturb), ::CompressibleEulerEquations2D) = ("rho", "v1", "v2", "theta")

###############################################################################
# Semidiscretization of the compressible Euler equations
density_current_setup = AtmoSetup()

equations = CompressibleEulerEquations2D(inertia_gravity_wave_setup.gamma)

boundary_conditions = (x_neg = boundary_condition_slip_wall,
                       x_pos = boundary_condition_slip_wall,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = FluxLMARS(340.0)

volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (25_600.0, 6_400.0)

cells_per_dimension = (128, 64)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (false, false))

semi = SemidiscretizationHyperbolic(mesh, equations, inertia_gravity_wave_setup, solver,
                                    source_terms = inertia_gravity_wave_setup,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 900.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "2_density_current/out",
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

visualization = VisualizationCallback(interval = 10000,
                                      solution_variables = cons2perturb,
                                      variable_names = variable_names)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback, visualization)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            maxiters = 1.0e7,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);

summary_callback()

# Visulization: 

# Compute the potential temperature perturbation
theta_perturb = let u = Trixi.wrap_array(sol.u[end], semi)
    @unpack g, c_p, c_v, p_0, theta_0 = density_current_setup
    rho = u[1, :, : ,:]
    rho_v1 = u[2, : ,: ,:]
    rho_v2 = u[3 , : ,: ,:]
    rho_e = u[4, : ,: ,:]
    R = c_p - c_v
    p = (equations.gamma - 1) .* (rho_e .- 0.5 .* (rho_v1 .* rho_v1 .+ rho_v2 .* rho_v2)./rho)
    rho_theta =  (p / p_0).^(c_v / c_p) .* p_0 / R
    rho_theta./rho .- theta_0
end

# Plotting the potential temperature perturbation
plot(ScalarPlotData2D(theta_perturb, semi), title = "Potential Temperature perturbation [K]")



