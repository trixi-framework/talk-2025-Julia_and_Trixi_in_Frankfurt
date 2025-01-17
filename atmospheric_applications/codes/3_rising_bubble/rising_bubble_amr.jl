using OrdinaryDiffEq
using Trixi
using Plots
using LaTeXStrings
# Definition of struct field to collect all the physical parameters and hydrostatic background state

struct WarmBubbleSetup
	# Physical constants
	g::Float64       # gravity of earth
	c_p::Float64     # heat capacity for constant pressure (dry air)
	c_v::Float64     # heat capacity for constant volume (dry air)
	gamma::Float64   # heat capacity ratio (dry air)
	p_0::Float64     # atmospheric pressure
	theta_0::Float64  # constant background potential temperature steady state     
	function WarmBubbleSetup(; g = 9.81, c_p = 1004.0, c_v = 717.0, gamma = c_p / c_v, p_0 = 100_000.0, theta_0 = 300.0)
		new(g, c_p, c_v, gamma, p_0, theta_0)
	end
end

# Initial condition for the rising bubble: constant potential temperature background state plus a perturbation
function (setup::WarmBubbleSetup)(x, t, equations::CompressibleEulerEquations2D)
	@unpack g, c_p, c_v, p_0, theta_0 = setup

	# center of perturbation
	center_x = 10000.0
	center_z = 2000.0
	# radius of perturbation
	radius = 2000.0
	# distance of current x to center of perturbation
	r = sqrt((x[1] - center_x)^2 + (x[2] - center_z)^2)

	# perturbation in potential temperature
	potential_temperature_ref = theta_0
	potential_temperature_perturbation = 0.0
	if r <= radius
		potential_temperature_perturbation = 2 * cospi(0.5 * r / radius)^2
	end
	potential_temperature = potential_temperature_ref + potential_temperature_perturbation

	# Exner pressure, solves hydrostatic equation for x[2]
	exner = 1 - g / (c_p * potential_temperature) * x[2]

	# pressure
	R = c_p - c_v    # gas constant (dry air)
	p = p_0 * exner^(c_p / R)

	# temperature
	T = potential_temperature * exner

	# density
	rho = p / (R * T)

	v1 = 20.0
	v2 = 0.0
	# E = c_v * T + 0.5 * (v1^2 + v2^2)
	# return SVector(rho, rho * v1, rho * v2, rho * E)
	return prim2cons(SVector(rho, v1, v2, p), equations)
end

# Source terms for the CompressibleEulerEquations2D
@inline function (setup::WarmBubbleSetup)(u, x, t, equations::CompressibleEulerEquations2D)
	@unpack g = setup
	rho, _, rho_v2, _ = u
	return SVector(zero(eltype(u)), zero(eltype(u)), -g * rho, -g * rho_v2)
end

@inline function cons2theta(u, equations::CompressibleEulerEquations2D)
	rho, rho_v1, rho_v2, rho_e = u
	@unpack g, c_p, c_v, p_0, theta_0 = warm_bubble_setup
	R = c_p - c_v
	p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * rho_v1 + rho_v2 * rho_v2) / rho)
	rho_theta = (p / p_0)^(c_v / c_p) * p_0 / R
	theta = rho_theta / rho - theta_0
	return theta
end

@inline function cons2perturb(u, equations::CompressibleEulerEquations2D)
	@unpack c_v, c_p, p_0, theta_0 = warm_bubble_setup
	R = c_p - c_v
	rho, rho_v1, rho_v2, rho_e = u
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
	rho_theta = (p / p_0)^(c_v / c_p) * p_0 / R
	theta = rho_theta / rho - theta_0

	return rho, v1, v2, theta

end

@inline Trixi.varnames(::typeof(cons2perturb), ::CompressibleEulerEquations2D) = ("rho", "v1", "v2", "theta")

###############################################################################
# Semidiscretization of the compressible Euler equations
warm_bubble_setup = WarmBubbleSetup()

equations = CompressibleEulerEquations2D(warm_bubble_setup.gamma)

boundary_conditions = Dict(
	:y_neg => boundary_condition_slip_wall,
	:y_pos => boundary_condition_slip_wall)


polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = FluxLMARS(340.0)

volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (20_000.0, 10_000.0)

trees_per_dimension = (4, 4)
mesh = P4estMesh(trees_per_dimension,
	polydeg = 3, initial_refinement_level = 2,
	coordinates_min = coordinates_min, coordinates_max = coordinates_max,
	periodicity = (true, false))
@show mesh
semi = SemidiscretizationHyperbolic(mesh, equations, warm_bubble_setup, solver,
	source_terms = warm_bubble_setup,
	boundary_conditions = boundary_conditions)
@show mesh
###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1000.0)  # 1000 seconds final time

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()


amr_indicator = IndicatorMax(semi, variable = cons2theta)
amr_controller = ControllerThreeLevel(semi, amr_indicator,
	base_level = 1,
	med_level = 3, med_threshold = 0.05,
	max_level = 4, max_threshold = 0.1)
amr_callback = AMRCallback(semi, amr_controller,
	interval = 5,
	adapt_initial_condition = true,
	adapt_initial_condition_only_refine = true)
    analysis_interval = 50
    alive_callback = AliveCallback(analysis_interval = analysis_interval)
    save_solution = SaveSolutionCallback(interval = analysis_interval,
    save_initial_solution = true,
    save_final_solution = true,
    output_directory = "3_rising_bubble/out",
    solution_variables = cons2perturb)



callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

###############################################################################
# run the simulation


##############
sol = solve(ode, SSPRK43(thread = OrdinaryDiffEq.True()),
	maxiters = 1.0e7,
	dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
	save_everystep = true, callback = callbacks);

summary_callback()

# Visulization: 

# Plotting the mesh
pd = PlotData2D(sol)
plot(getmesh(pd))


# Compute the potential temperature perturbation
theta_perturb = let u = Trixi.wrap_array(sol.u[end], semi)
	@unpack g, c_p, c_v, p_0, theta_0 = warm_bubble_setup
	rho = u[1, :, :, :]
	rho_v1 = u[2, :, :, :]
	rho_v2 = u[3, :, :, :]
	rho_e = u[4, :, :, :]
	R = c_p - c_v
	p = (equations.gamma - 1) .* (rho_e .- 0.5 .* (rho_v1 .* rho_v1 .+ rho_v2 .* rho_v2) ./ rho)
	rho_theta = (p / p_0) .^ (c_v / c_p) .* p_0 / R
	rho_theta ./ rho .- theta_0
end

# Plotting the potential temperature perturbation
plot(ScalarPlotData2D(theta_perturb, semi), title = "Potential Temperature perturbation [K]")