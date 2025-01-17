### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ fe03cfcd-bd4b-4c77-99b4-719f51dbe893
begin
	using Trixi
	using OrdinaryDiffEq
	using Plots
	using LaTeXStrings
end

# ╔═╡ ad6ddb85-006a-4ee6-92da-4f56ea4780b0
md""" 
# Atmospheric Applications in Trixi.jl

Some famous benchmark test cases will be illustrated, showing different possibilities and capabalities of Trixi.jl.

- Setting up your own problem: Rising Bubble
- Warped/Curvilinear Mesh
- Adaptive Mesh Refinement (AMR)
- Mountain Waves: Orography and HOHQ Mesh
"""

# ╔═╡ 065746f1-1deb-49c4-854d-457a8f7919c2
md""" 
# Setting up your own problem in Trixi.jl

Simulations in Trixi.jl consists of different building blocks
- Set of equations
- Source terms
- Initial conditions
- Boundary Conditions
- Domain and Mesh (P4est, T8code, StructuredMesh ...)
- Semidiscretization (DGSEM, FDSBP, ...)
- Time Integrator

"""

# ╔═╡ 2f015b91-8c8a-42c7-8264-d4c722591b64
md""" 
# Rising Bubble test case

In Trixi.jl you can choose among a broad set of equations. For dry air atmospheric applications we choose 2D Compressible Euler Equations.
"""

# ╔═╡ 482f659c-ffbe-44db-bf95-a7bdd86c0640
gamma = 1004.0/717.0

# ╔═╡ e1947159-3084-472c-b9a5-fbdf648042ef
equations  = CompressibleEulerEquations2D(gamma)

# ╔═╡ 5c644ba6-f209-4336-92aa-170164bbe99c
md""" 
#### Gravity source term
Pointwise source terms can be easily added to the right-hand side of the equations with user-defined functions or choosing among the built-in source terms functions in Trixi.jl
"""

# ╔═╡ 287e2d6c-3b9c-4ba9-9d01-42f014047d81
md"
```math
\begin{equation}
    \partial_t \begin{pmatrix}
           \varrho \\
           \varrho u \\
          \varrho v \\
           \varrho E
         \end{pmatrix} + \partial_x \begin{pmatrix}
           \varrho u \\
           \varrho u^2 + p \\
           \varrho u v \\
            (\varrho E + p) u
         \end{pmatrix}+ \partial_y \begin{pmatrix}
           \varrho u \\
           \varrho u v \\
           \varrho v^2 + p\\
            (\varrho E + p) v
         \end{pmatrix} =  \begin{pmatrix}
           0 \\
           0 \\
           - \varrho g \\
           -\varrho g v 
         \end{pmatrix}
\end{equation}
```
"

# ╔═╡ 012b2045-4172-41f6-98bb-c773b92cf1c3
function source_terms_gravity(u, x, t, equations::CompressibleEulerEquations2D)
    g = 9.81
	rho, _, rho_v2, _ = u
    return SVector(0, 0, -g * rho, -g * rho_v2)
end

# ╔═╡ f94051e4-a809-4732-9afe-550c3fb43e44
md""" 
#### Initial condition: warm rising bubble
In an hydrostatic balance steady state condition with constant potential temperature $\theta = 300$, the motion is driven by the following perturbation.
"""

# ╔═╡ 3fa067c0-7491-4a55-8b9b-16376a48a90b
L"""
\theta' =
\begin{cases}
0, & \text{se } r > r_c, \\
\frac{\theta_c}{2} \left[ 1 + \cos\left(\frac{\pi r}{r_c}\right) \right], & \text{se } r \leq r_c.
\end{cases}
"""

# ╔═╡ 61de0650-e3cb-4dfe-823b-65ca76b843dd
function initial_condition_rising_bubble(x, t, equations::CompressibleEulerEquations2D)
	g = 9.81
	c_p = 1004.0
	c_v = 717.0
	p_0 = 100_000.0
	R = c_p - c_v
	potential_temperature_ref = 300.0
	# center of perturbation
    center_x = 10000.0
    center_z = 2000.0
    # radius of perturbation
    radius = 2000.0
    # distance of current x to center of perturbation
    r = sqrt((x[1] - center_x)^2 + (x[2] - center_z)^2)

    # perturbation in potential temperature
    potential_temperature_perturbation = 0.0
    if r <= radius
        potential_temperature_perturbation = 2 * cospi(0.5 * r / radius)^2
    end
    potential_temperature = potential_temperature_ref + 	   potential_temperature_perturbation

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 - g / (c_p * potential_temperature) * x[2]

    # pressure
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

# ╔═╡ ea43bccf-c1c0-4661-8438-6bd2466d411b
md""" 
#### Boundary conditions
Along the $x$-direction periodic boundary conditions are employed and for the top and bottom walls the boundary conditions are no-flux.
"""

# ╔═╡ 34f5fbd0-abe8-468d-9a62-30d5e5924b97
boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall); nothing

# ╔═╡ f583fefc-d0a6-4586-8c48-616942d29e98
md""" 
#### Defining a simple Mesh

The mesh is a rectangular domain of size $[0, 20 \text{ km}] \text{ x } [0, 10 \text{ km}]$.

The domain is divided into 64 cells along the $x$ horizontal direction and 32 in the $z$ vertical direction.
"""

# ╔═╡ 022ee8fb-c758-4f87-8b60-8103cd8719f5
coordinates_min = (0.0, 0.0)

# ╔═╡ b1290b6e-8a81-41c9-9249-b00d5bfbfd5d
coordinates_max = (20_000.0, 10_000.0)

# ╔═╡ 4e055027-c7b4-4086-8a8c-a19bf97769c5
cells_per_dimension = (64, 32)

# ╔═╡ 0e50dafd-a179-4c22-8afa-fa96cd45ff43
mesh = StructuredMesh(cells_per_dimension, 
					  coordinates_min, 
					  coordinates_max,
                      periodicity = (true, false)); nothing

# ╔═╡ d53556c2-a1d3-48c4-a265-2d0a61b5ec9a
md""" 
#### DGSEM Discretization

DGSEM Flux differencing semi-discretization has been chosen for these test cases. 
- Third order polynomial degree
- LMARS flux for the surface integral
- Kennedy-Gruber flux for the volume integral (symmetric!)
"""

# ╔═╡ ff35fbf5-53fd-4e9f-9090-5aab89812b31
begin
	polydeg = 3
	
	surface_flux = FluxLMARS(340.0)
	
	volume_flux = flux_kennedy_gruber
	
	volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
	
	solver = DGSEM(polydeg, surface_flux, volume_integral); nothing
end

# ╔═╡ 00f13c5c-9c81-4a68-8885-473d3d61bbf4
semi = SemidiscretizationHyperbolic(mesh, 
									equations, 
									initial_condition_rising_bubble, 
									solver,
                                    source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions); nothing

# ╔═╡ fdd7e368-c826-49ac-aa9a-5002d0a0cae7
begin
function plot_mesh(semi)	

tspan = (0.0, 0.0)
ode = semidiscretize(semi, tspan)
sol = solve(ode, SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false);	
pd = PlotData2D(sol)
plot(getmesh(pd))
end
	plot_mesh(semi)
end

# ╔═╡ 22a5f545-e62e-4594-92ba-7d7c90d172d4
md""" 
#### Creating ODE function and defining the time integrator.

Once the ODE has been built, we are ready to pass our RHS function to the ODE solver. The simulation is running with SSPRK43 with multiple-threads for $T = 1000$ s physical time.
"""

# ╔═╡ f5691db4-fbba-4356-8898-63a33fe56f10
tspan = (0.0, 1000.0)  # 1000 seconds final time

# ╔═╡ bbdb5d15-b207-4b62-a2d6-d496ec8af214
ode = semidiscretize(semi, tspan); nothing

# ╔═╡ d2db4dde-821c-4a4c-80d3-a172b777782a
sol = solve(ode, 
		    SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false);

# ╔═╡ f4a2de9f-357c-424f-9619-500d3e524c3f
# ╠═╡ disabled = true
#=╠═╡
begin
pd = PlotData2D(sol)
plot(getmesh(pd))
end
  ╠═╡ =#

# ╔═╡ 14670101-6792-4587-8efb-865627b4e446
md"### Initial condition rising bubble"

# ╔═╡ 54b02548-63e0-4b07-b982-56c51e6f450a
begin
theta_perturb_0 = let u = Trixi.wrap_array(sol.u[1], semi)
    rho = u[1, :, : ,:]
    rho_v1 = u[2, : ,: ,:]
    rho_v2 = u[3 , : ,: ,:]
    rho_e = u[4, : ,: ,:]
	c_p = 1004.0
	c_v = 717.0
    R = c_p - c_v
	p_0 = 100_000.0
	potential_temperature_ref = 300.0
    p = (equations.gamma - 1) .* (rho_e .- 0.5 .* (rho_v1 .* rho_v1 .+ rho_v2 .* rho_v2)./rho)
    rho_theta =  (p / p_0).^(c_v / c_p) .* p_0 / R
    rho_theta./rho .- potential_temperature_ref
end

plot(ScalarPlotData2D(theta_perturb_0, semi), title = "Potential Temperature perturbation [K], t = 0")
end

# ╔═╡ 695bf5b1-2f99-4f92-8c35-2eb0279a5663
md"### Solution at t = 1000s  rising bubble"

# ╔═╡ e6cdc448-565e-4e50-a726-75847ef13bc5
begin
theta_perturb = let u = Trixi.wrap_array(sol.u[end], semi)
    rho = u[1, :, : ,:]
    rho_v1 = u[2, : ,: ,:]
    rho_v2 = u[3 , : ,: ,:]
    rho_e = u[4, : ,: ,:]
	c_p = 1004.0
	c_v = 717.0
    R = c_p - c_v
	p_0 = 100_000.0
	potential_temperature_ref = 300.0
    p = (equations.gamma - 1) .* (rho_e .- 0.5 .* (rho_v1 .* rho_v1 .+ rho_v2 .* rho_v2)./rho)
    rho_theta =  (p / p_0).^(c_v / c_p) .* p_0 / R
    rho_theta./rho .- potential_temperature_ref
end

plot(ScalarPlotData2D(theta_perturb, semi), title = "Potential Temperature perturbation [K], t = 1000 s")
end

# ╔═╡ bb0b8663-98a9-44ab-ac7c-4c10ba56378e
md"### Elixir Summary"

# ╔═╡ 83825531-2488-4133-89c8-a88743cf10b1
begin 
	function simple_elixir()
      # Equations
	  equations  = CompressibleEulerEquations2D(gamma)

	  # Boundary conditions	 
	  boundary_conditions = (x_neg = boundary_condition_periodic,
                             x_pos = boundary_condition_periodic,
                             y_neg = boundary_condition_slip_wall,
                             y_pos = boundary_condition_slip_wall)

      # Initial condition	
	  initial_condition = initial_condition_rising_bubble	

      # Source terms	
	  source_terms = source_terms_gravity
	
	  # Mesh
	  coordinates_min = (0.0, 0.0)
      coordinates_max = (20_000.0, 10_000.0)

	  cells_per_dimension = (32, 32)
      mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))	
	
	  # DGSEM semidiscretization
	  polydeg = 3
	
	  surface_flux = FluxLMARS(340.0)
	
	  volume_flux = flux_kennedy_gruber
	
	  volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
	
	  solver = DGSEM(polydeg, surface_flux, volume_integral)
	  
	  semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)
	  # ODE Solver
	  tspan = (0.0, 1000.0)
	
      ode = semidiscretize(semi, tspan)

	  sol = solve(ode, 
		    SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false)
	  end
end

# ╔═╡ 6a6fd640-396e-41fa-8c6b-977144181ad1
md""" 
## Warped-Curvilinear Mesh.

We consider the same mesh and apply a warping transformation.

"""

# ╔═╡ 8fc7fe79-5b78-49a9-9bca-2179460324cc
function mapping(xi, eta)
    x = xi + 0.1 * sin(2 * pi * xi) * sin(2 * pi * eta)
    y = eta + 0.1 * sin(2 * pi * xi) * sin(2 * pi * eta)
    return SVector(20_000.0*0.5 * (1 + x), 10_000.0*0.5 * (1+y))
end

# ╔═╡ 5bb04fa4-36f4-4207-82f1-33a64797ceae
mesh_warped = StructuredMesh(cells_per_dimension, mapping,
                             periodicity = (true, false)); nothing

# ╔═╡ 956d6c3d-0138-4791-95a9-9f8525570a32
begin
semi_warped = SemidiscretizationHyperbolic(mesh_warped, 
									equations, 
									initial_condition_rising_bubble, 
									solver,
                                    source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions)
tspan_warped = (0.0, 0.0)
ode_warped = semidiscretize(semi_warped, tspan_warped)
sol_warped = solve(ode_warped, SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false);	
pd_warped = PlotData2D(sol_warped)
plot(getmesh(pd_warped))
end

# ╔═╡ 0fc2fdc1-4d45-428f-9a47-590409646921
md""" 
## Adaptive-Mesh Refinement

Simulations in Trixi.jl can be sped up maintaining an overall good accuracy through Adaptive Mesh Refinement. The hierarchical Cartesian mesh is locally refined, based on a chosen reference variable.

In Trixi.jl some callback functions can be defined, for numerical analysis and step size control. The AMR function is also defined passing an AMR Callback to the ODE solver, as shown in this example.

For AMR we switch from StructuredMesh to P4estMesh.
"""

# ╔═╡ 98a46c89-8805-4461-8528-3bd5f4433e31
begin
	trees_per_dimension = (4, 4)
	mesh_amr = P4estMesh(trees_per_dimension,
						 polydeg = 3, initial_refinement_level = 2,
		                 coordinates_min = coordinates_min, 
		                 coordinates_max = coordinates_max,
		                 periodicity = (true, false)); nothing
end

# ╔═╡ d9f92197-b9dc-4820-aca0-ad3dd2e30080
md"""

The semi-discretization follows the same step as before.



"""

# ╔═╡ 9051a970-0f5f-4019-b8ee-37b44753865b
begin
boundary_conditions_amr = Dict( :y_neg => boundary_condition_slip_wall,
							    :y_pos => boundary_condition_slip_wall )

semi_amr = SemidiscretizationHyperbolic(mesh_amr, 
										equations, 
									 	initial_condition_rising_bubble, 
										solver,
                                    	source_terms = source_terms_gravity,
                                    	boundary_conditions = 			           boundary_conditions_amr)
ode_amr = semidiscretize(semi_amr, (0.0, 1000.0)); nothing	
end

# ╔═╡ 498d3262-df89-408e-8ede-f32268894709
function cons2theta(u, equations::CompressibleEulerEquations2D)
	rho, rho_v1, rho_v2, rho_e = u
	g = 9.81
	c_p = 1004.0
	c_v = 717.0
	p_0 = 100_000.0
	theta0 = 300.0
	R = c_p - c_v
	p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * rho_v1 + rho_v2 * rho_v2) / rho)
	rho_theta = (p / p_0)^(c_v / c_p) * p_0 / R
	theta = rho_theta / rho - theta0
	return theta
end ;nothing

# ╔═╡ 28175473-1386-483f-9443-817cad9fd53e
begin
amr_indicator = IndicatorMax(semi_amr, variable = cons2theta)

amr_controller = ControllerThreeLevel(semi_amr, amr_indicator,
	base_level = 1,
	med_level = 3, med_threshold = 0.05,
	max_level = 4, max_threshold = 0.1)
	
amr_callback = AMRCallback(semi_amr, amr_controller,
	interval = 5,
	adapt_initial_condition = true,
	adapt_initial_condition_only_refine = true)
end

# ╔═╡ 7ea5ea01-2cbc-4f0c-9d8b-64a1fc83294c
md"""Once the callback is set, it can be passed to the ODE solver as a key-word argument. The mesh and solution can be plotted to see where the mesh has been refined. On my machine the solution takes almost half the time for the refined mesh and it shows visibly better results than the static mesh."""

# ╔═╡ f2d35fb1-cba0-419b-b3fc-3fee9a152532
callbacks_amr = CallbackSet(amr_callback)

# ╔═╡ 029a886b-edec-4238-85a2-39f10dde3f5f
sol_amr = solve(ode_amr, SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false, callback = callbacks_amr);

# ╔═╡ 7ebed9ca-8f57-4971-a274-bb816437642f
begin
pd_amr = PlotData2D(sol_amr)

theta_perturb_amr = let u = Trixi.wrap_array(sol_amr.u[end], semi_amr)
    rho = u[1, :, : ,:]
    rho_v1 = u[2, : ,: ,:]
    rho_v2 = u[3 , : ,: ,:]
    rho_e = u[4, : ,: ,:]
	c_p = 1004.0
	c_v = 717.0
    R = c_p - c_v
	p_0 = 100_000.0
	potential_temperature_ref = 300.0
    p = (equations.gamma - 1) .* (rho_e .- 0.5 .* (rho_v1 .* rho_v1 .+ rho_v2 .* rho_v2)./rho)
    rho_theta =  (p / p_0).^(c_v / c_p) .* p_0 / R
    rho_theta./rho .- potential_temperature_ref
end

plot(ScalarPlotData2D(theta_perturb_amr, semi_amr), title = "Potential Temperature perturbation [K], t = 1000 s")

plot!(getmesh(pd_amr))
end

# ╔═╡ c886e8ca-8506-4848-a342-2613747cddb1
begin
	plot(ScalarPlotData2D(theta_perturb_amr, semi_amr), title = "Potential Temperature perturbation [K], t = 1000 s")
end

# ╔═╡ ef5bba7b-2406-4e82-b509-8cace58bf508
md"""

# Hydrostatic and non-Hydrostatic waves

The versiera di Agnesi mountain profile is used for the following test case.

$h(x, z) = \frac{h_c}{1 + \left( \frac{x-x_c}{a_c}\right)^2}$

The solution evolves until a steady-state solution of a linear hydrostatic flow is reached over a single-peaked mountain. The domain size is $[-120 \text{ km}, 120 \text{ km}] \text{ x } [0, 30 \text{ km}]$. 

The domain is divided into 32 cells along the $x$ horizontal direction and 32 in the $z$ vertical direction.

Below a representation of the versiera di Agnesi mountain profile mesh with a peak of 5km.
"""

# ╔═╡ d74493e3-747a-42ac-86ec-3f2b0d01eabb
begin
	function mesh_agnesi_profile(; peak = 1.0)
	a = 10000.0
	L = 240000.0
	H = 30000.0
	y_b = peak / (1 + (L/2 / a)^2)
	alfa = (H - y_b) * 0.5

	f1(s) = SVector(-L/2, y_b + alfa * (s + 1)) 				# left
	f2(s) = SVector(L/2, y_b + alfa * (s + 1))  				# right
	f3(s) = SVector(s * L/2, peak / (1 + ( s * L/2)^2 / a^2))	# bottom
	f4(s) = SVector(s * L/2, H) 								# top

	mesh = StructuredMesh((32, 32), (f1, f2, f3, f4), periodicity = (true, false))
		return mesh
	end
	mesh_orography_v = mesh_agnesi_profile(peak = 5000.0); nothing
end

# ╔═╡ fc94d662-30cd-4ef9-9a2e-3df9c4e8a498
md""" 
#### Rayleigh damping profiles

On the right-hand side to avoid reflective waves on the wall that introduces noise in the numerical simulation and makes it spurious, an absorbing sponge layer which relaxes the numerical solution is prescribed. The waves are damped in the last 10 km of the top layer.

$S(z) = - \frac{\alpha}{2} \left ( 1 - \cos\left (\pi \frac{z - z_B}{z_T - z_B}\right ) \right)$

"""

# ╔═╡ 5b9eaae3-9278-485d-a119-40004f7bfaec
function source_terms_damping(u, x, t, equations::CompressibleEulerEquations2D)
	c_p = 1004.0; c_v = 717.0; p_0 = 100_000.0
	gamma = equations.gamma
	u0 = 20.0
	z_B = 15000.0; z_T = 30000.0
	T_0 = 250.0
	g = 9.81
	R = c_p - c_v
	Nf = g/sqrt(c_p*T_0)
	rho, rho_v1, rho_v2, rho_e = u
    
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
	rho_theta =  (p / p_0)^(c_v / c_p) * p_0 / R
    theta = rho_theta/rho

    alfa = 0.1

    if x[2] <= z_B
        S_v = 0.0
    elseif (x[2] - z_B)/(z_T - z_B) <= 1/2
        S_v = -alfa/2 *(1 - cospi((x[2] - z_B)/(z_T - z_B)))
    else 
        S_v = -alfa/2 * ( 1 + ((x[2] - z_B)/(z_T - z_B) - 1/2)) * pi
    end

    exner = exp(-Nf^2/g * x[2])

    theta_0 = T_0/exner
    K = p_0 * (R/p_0)^gamma
	du2 = rho * (v1-u0) * S_v
	du3 = rho_v2 * S_v 
	du4 = rho * (theta-theta_0) * S_v * K * gamma/(gamma - 1.0) * (rho_theta)^(gamma - 1.0)  + du2 * v1 + du3 * v2

	return SVector(zero(eltype(u)), du2, du3 -g * rho, du4 - g * rho_v2)

end

# ╔═╡ ebecb24d-502e-4358-828f-7b3d03129d5f
function initial_condition_linear_hydrostatic(x, t, equations::CompressibleEulerEquations2D)
	g = 9.81
	c_p = 1004.0
	c_v = 717.0
	p_0 = 100_000.0
	T_0 = 250.0
	u0 = 20.0
	Nf = g/sqrt(c_p*T_0)

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = exp(-Nf^2/g*x[2])
    # pressure
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)

    # density
    rho = p / (R * T_0)
    v1 = u0
    v2 = 0.0
    
    return prim2cons(SVector(rho, v1, v2 ,p), equations)
end; nothing

# ╔═╡ e507bb22-287a-46ac-9932-745a103578d2
begin
semi_hydrostatic_v = SemidiscretizationHyperbolic(mesh_orography_v, 
												equations, initial_condition_linear_hydrostatic,solver, 
												source_terms = source_terms_damping,
                                                boundary_conditions = 			      												boundary_conditions)
tspan_hydrostatic_v = (0.0, 0.0)
ode_hydrostatic_v = semidiscretize(semi_hydrostatic_v, tspan_hydrostatic_v)
sol_hydrostatic_v = solve(ode_hydrostatic_v, SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false);	
pd_hydrostatic_v = PlotData2D(sol_hydrostatic_v)
plot(getmesh(pd_hydrostatic_v), aspect_ratio = 10)
end

# ╔═╡ 2b86d029-f10d-46c4-a1e9-28fbb234663d
mesh_orography = mesh_agnesi_profile(); nothing

# ╔═╡ 85431810-87ee-4dec-bbfa-dfff8c5eda04
begin
semi_hydrostatic = SemidiscretizationHyperbolic(mesh_orography, 
												equations, initial_condition_linear_hydrostatic,solver, 
												source_terms = source_terms_damping,
                                                boundary_conditions = 			      												boundary_conditions)
tspan_hydrostatic = (0.0, 5*3600.0)
ode_hydrostatic = semidiscretize(semi_hydrostatic, tspan_hydrostatic)
sol_hydrostatic = solve(ode_hydrostatic, SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false);	
pd_hydrostatic = PlotData2D(sol_hydrostatic); nothing
end

# ╔═╡ 492a2130-8bf5-4a00-adf5-8f39656cebda
md""" ###### Solution linear hydrostatic mountain at t = 5h"""

# ╔═╡ 9455e131-3968-4d57-a9f9-3b68bd9dcac8
begin
	horizontal_velocity = let u = Trixi.wrap_array(sol_hydrostatic.u[end], semi_hydrostatic)
	    rho = u[1, :, : ,:]
	    rho_v1 = u[2, : ,: ,:]
	    rho_v1 ./ rho .- 20.0
	end
	# Convertire i limiti in km
xtick_positions = [-40000, -20000, 20000, 40000]  # Posizioni dei tick in metri
ytick_positions = [0, 2000, 4000, 6000, 8000, 10000, 12000]          # Posizioni dei tick in metri

xtick_labels = ["-40 km", "-20 km", "20 km", "40 km"]  # Etichette da visualizzare
ytick_labels = ["0 km", "2 km", "4 km", "6 km", "8 km", "10 km", "12 km"]       
	# Plotting the vertical velocity component
	plot(ScalarPlotData2D(horizontal_velocity, semi_hydrostatic), title = "Horizontal velocity component [m/s]", aspect_ratio = 5, xlim = (-40000, 40000), ylim = (0, 12000), 
     xticks = (xtick_positions, xtick_labels),  # Specifica i tick x
     yticks = (ytick_positions, ytick_labels))
end

# ╔═╡ 2617f478-5556-4ab1-a506-a0a7bc0238cc
md""" 
## Mountain Schär

The versiera di Agnesi mountain profile is used for the following test case.

$h(x, z) = h_c e^{-\frac{x^2}{a_c^2}} \cos^2\left(\frac{\pi x}{\lambda_c} \right )$

The solution evolves until a steady-state solution of a linear hydrostatic flow is reached, over a single-peaked mountain. The domain size is $[-25, \text{ }25 ] \text{ x } [0, \text{ }21] \text{ km}$. 
"""

# ╔═╡ e44793c5-aa7c-4b64-bdc7-52a6bd60980c
begin
function mountain_schar_profile()
	a = 5000.0
	L = 50000.0
	H = 30000.0
	lambda_c = 4000.0
	hc = 250.0
	y_b = hc * exp(-(L/2/a)^2)*cospi(L/2/lambda_c)^2
	alfa = (H - y_b) * 0.5

	f1(s) = SVector(-L/2, y_b + alfa * (s + 1))
	f2(s) = SVector(L/2, y_b + alfa * (s + 1))
	f3(s) = SVector(s * L/2, hc * exp(-(s * L/2 /a)^2) * cospi(s * L/2 /lambda_c)^2)
	f4(s) = SVector(s * L/2, H)

	mesh = StructuredMesh((32, 32), (f1, f2, f3, f4), periodicity = (true, false))
	
	return mesh
end
	mesh_schar = mountain_schar_profile(); nothing
end

# ╔═╡ a45a1305-faa2-4c6e-a79b-9c771e7ff7ec
function source_terms_rayleigh(u, x, t, equations::CompressibleEulerEquations2D)
	g = 9.81
	c_p = 1004.0
	c_v = 717.0
	gamma = c_p/c_v
	p_0 = 100_000.0
	theta_0 = 280.0
	z_B = 15000.0
	z_T = 21000.0
	Nf = 0.01
	u0 = 10.0
	R = c_p - c_v
	
	rho, rho_v1, rho_v2, rho_e = u
    
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
	rho_theta =  (p / p_0)^(c_v / c_p) * p_0 / R
    theta = rho_theta/rho

    alfa = 0.1

    x_B = 20000.0
    x_T = 25000.0

    if x[2] <= z_B
        S_v = 0.0
    else
        S_v = -alfa * sinpi(0.5*(x[2] - z_B)/(z_T - z_B))^2
    end
    if x[1] < x_B
        S_h1 = 0.0
    else
        S_h1 = -alfa * sinpi(0.5*(x[1] - x_B)/(x_T - x_B))^2
    end

    if x[1] > -x_B
        S_h2 = 0.0
    else
        S_h2 = -alfa * sinpi(0.5*(x[1] + x_B)/(-x_T + x_B))^2
    end

    K = p_0 * (R/p_0)^gamma
	du2 = rho * (v1-u0) * (S_v + S_h1 + S_h2)
	du3 = rho_v2 * (S_v + S_h1 + S_h2) 
	du4 = rho * (theta-theta_0) * (S_v + S_h1 + S_h2) * K * gamma/(gamma - 1.0) * (rho_theta)^(gamma - 1.0)  + du2 * v1 + du3 * v2

	return SVector(zero(eltype(u)), du2, du3 - g * rho, du4*0.0 - g * rho_v2)

end; nothing

# ╔═╡ 19818181-c2f4-4824-8eed-f938a3b204cf
function initial_condition_schar(x, t, equations::CompressibleEulerEquations2D)
	g = 9.81; c_p = 1004.0; c_v = 717.0; p_0 = 100_000.0; theta_0 = 280.0; u0 = 10.0
	Nf = 0.01
    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 + g^2 / (c_p * theta_0 * Nf^2) * (exp(-Nf^2 / g * x[2]) - 1)
    # pressure
    p_0 = 100_000.0  # reference pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p_0 * exner^(c_p / R)
    potential_temperature = theta_0 * exp(Nf^2/g*x[2])
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = u0
    v2 = 0.0
    
    return prim2cons(SVector(rho, v1, v2 ,p), equations)
end; nothing

# ╔═╡ e7958ec8-b910-4e03-90bf-f5e50d5af4e7
begin
	semi_schar = SemidiscretizationHyperbolic(mesh_schar, equations, initial_condition_schar, solver, source_terms = source_terms_rayleigh,
                                    boundary_conditions = boundary_conditions)
	
	tspan_schar = (0.0, 1*3600.0)
	ode_schar = semidiscretize(semi_schar, tspan_schar)
	sol_schar = solve(ode_schar, 
		              SSPRK43(thread = OrdinaryDiffEq.True()),
                      maxiters = 1.0e7,
                      dt = 1.0,
                      save_everystep = false); nothing	
end

# ╔═╡ f91c4f28-38b1-4267-8c97-d494f4535cb9
begin
	u1 = let u = Trixi.wrap_array(sol_schar.u[end], semi_schar)
	    rho = u[1, :, : ,:]
	    rho_v1 = u[2, : ,: ,:]
	    rho_v1 ./ rho .- 20.0
	end
	function plot_schar(u1, semi_schar)
	xtick_positions = [-10000, -5000, 5000, 10000]  
	ytick_positions = [0, 2000, 4000, 6000, 8000, 10000]

	xtick_labels = ["-10 km", "-5 km", "5 km", "10 km"] 
	ytick_labels = ["0 km", "2 km", "4 km", "6 km", "8 km", "10 km"]       
	
	plot(ScalarPlotData2D(u1, semi_schar), title = "Horizontal velocity component [m/s]", aspect_ratio = 2, xlim = (-10000, 10000),  ylim = (0, 10000), 
     xticks = (xtick_positions, xtick_labels), 
     yticks = (ytick_positions, ytick_labels))
	end
	plot_schar(u1, semi_schar)
end

# ╔═╡ 77118cea-fd54-4954-944b-c13ed6e0ad1d
pd_schar = PlotData2D(sol_schar); nothing

# ╔═╡ 5533f9d4-8ce1-489b-adad-6016bc003e2e
begin
	function plot_mesh_schar(pd_schar)
		# Convertire i limiti in km
		xtick_positions = [-10000, -5000, 5000, 10000]  # Posizioni dei tick in metri
		ytick_positions = [0, 2000, 4000, 6000, 8000, 10000]          # Posizioni dei tick in metri
	
		xtick_labels = ["-10 km", "-5 km", "5 km", "10 km"]  # Etichette da visualizzare
		ytick_labels = ["0 km", "2 km", "4 km", "6 km", "8 km", "10 km"]       
		# Plotting the vertical velocity component
		plot(getmesh(pd_schar), aspect_ratio = 2, xlim = (-10000, 10000),  ylim = (0, 10000), 
	     xticks = (xtick_positions, xtick_labels),  # Specifica i tick x
	     yticks = (ytick_positions, ytick_labels))
	end
	
	plot_mesh_schar(pd_schar)
end

# ╔═╡ 8be590dc-9b54-4454-88fc-1fa82e91cd5d
md"""
# HOHQ Mesh

[HOHQMesh.jl](https://github.com/trixi-framework/HOHQMesh.jl)


HOHQ Mesh is a library that allows the user to define high order hex/quadrilater meshes even for complex geometries that cannot be easily described by analytical functions, such as irregular surfaces.

As an example, here is shown a generic mountain profile given by single data point and interpolated. The code is available in the repository.
"""

# ╔═╡ 3427fa7c-41f6-4b08-af59-3eb9d72a7f1a
md"""
# Conclusions

- Implementing your own set of equations.

- SVG to mesh via HOHQ Mesh

- [TrixiAtmo.jl](https://github.com/trixi-framework/TrixiAtmo.jl)

"""

# ╔═╡ Cell order:
# ╟─ad6ddb85-006a-4ee6-92da-4f56ea4780b0
# ╠═fe03cfcd-bd4b-4c77-99b4-719f51dbe893
# ╟─065746f1-1deb-49c4-854d-457a8f7919c2
# ╟─2f015b91-8c8a-42c7-8264-d4c722591b64
# ╟─482f659c-ffbe-44db-bf95-a7bdd86c0640
# ╠═e1947159-3084-472c-b9a5-fbdf648042ef
# ╟─5c644ba6-f209-4336-92aa-170164bbe99c
# ╟─287e2d6c-3b9c-4ba9-9d01-42f014047d81
# ╠═012b2045-4172-41f6-98bb-c773b92cf1c3
# ╟─f94051e4-a809-4732-9afe-550c3fb43e44
# ╟─3fa067c0-7491-4a55-8b9b-16376a48a90b
# ╠═61de0650-e3cb-4dfe-823b-65ca76b843dd
# ╟─ea43bccf-c1c0-4661-8438-6bd2466d411b
# ╠═34f5fbd0-abe8-468d-9a62-30d5e5924b97
# ╠═f583fefc-d0a6-4586-8c48-616942d29e98
# ╟─022ee8fb-c758-4f87-8b60-8103cd8719f5
# ╟─b1290b6e-8a81-41c9-9249-b00d5bfbfd5d
# ╟─4e055027-c7b4-4086-8a8c-a19bf97769c5
# ╠═0e50dafd-a179-4c22-8afa-fa96cd45ff43
# ╟─fdd7e368-c826-49ac-aa9a-5002d0a0cae7
# ╟─d53556c2-a1d3-48c4-a265-2d0a61b5ec9a
# ╠═ff35fbf5-53fd-4e9f-9090-5aab89812b31
# ╠═00f13c5c-9c81-4a68-8885-473d3d61bbf4
# ╟─22a5f545-e62e-4594-92ba-7d7c90d172d4
# ╟─f5691db4-fbba-4356-8898-63a33fe56f10
# ╠═bbdb5d15-b207-4b62-a2d6-d496ec8af214
# ╠═d2db4dde-821c-4a4c-80d3-a172b777782a
# ╟─f4a2de9f-357c-424f-9619-500d3e524c3f
# ╟─14670101-6792-4587-8efb-865627b4e446
# ╟─54b02548-63e0-4b07-b982-56c51e6f450a
# ╟─695bf5b1-2f99-4f92-8c35-2eb0279a5663
# ╟─e6cdc448-565e-4e50-a726-75847ef13bc5
# ╟─bb0b8663-98a9-44ab-ac7c-4c10ba56378e
# ╠═83825531-2488-4133-89c8-a88743cf10b1
# ╟─6a6fd640-396e-41fa-8c6b-977144181ad1
# ╠═8fc7fe79-5b78-49a9-9bca-2179460324cc
# ╠═5bb04fa4-36f4-4207-82f1-33a64797ceae
# ╟─956d6c3d-0138-4791-95a9-9f8525570a32
# ╠═0fc2fdc1-4d45-428f-9a47-590409646921
# ╠═98a46c89-8805-4461-8528-3bd5f4433e31
# ╟─d9f92197-b9dc-4820-aca0-ad3dd2e30080
# ╠═9051a970-0f5f-4019-b8ee-37b44753865b
# ╠═28175473-1386-483f-9443-817cad9fd53e
# ╟─498d3262-df89-408e-8ede-f32268894709
# ╟─7ea5ea01-2cbc-4f0c-9d8b-64a1fc83294c
# ╠═f2d35fb1-cba0-419b-b3fc-3fee9a152532
# ╠═029a886b-edec-4238-85a2-39f10dde3f5f
# ╟─7ebed9ca-8f57-4971-a274-bb816437642f
# ╟─c886e8ca-8506-4848-a342-2613747cddb1
# ╠═ef5bba7b-2406-4e82-b509-8cace58bf508
# ╟─d74493e3-747a-42ac-86ec-3f2b0d01eabb
# ╟─e507bb22-287a-46ac-9932-745a103578d2
# ╟─fc94d662-30cd-4ef9-9a2e-3df9c4e8a498
# ╠═5b9eaae3-9278-485d-a119-40004f7bfaec
# ╟─ebecb24d-502e-4358-828f-7b3d03129d5f
# ╟─2b86d029-f10d-46c4-a1e9-28fbb234663d
# ╟─85431810-87ee-4dec-bbfa-dfff8c5eda04
# ╠═492a2130-8bf5-4a00-adf5-8f39656cebda
# ╟─9455e131-3968-4d57-a9f9-3b68bd9dcac8
# ╠═2617f478-5556-4ab1-a506-a0a7bc0238cc
# ╠═e44793c5-aa7c-4b64-bdc7-52a6bd60980c
# ╟─5533f9d4-8ce1-489b-adad-6016bc003e2e
# ╟─a45a1305-faa2-4c6e-a79b-9c771e7ff7ec
# ╟─19818181-c2f4-4824-8eed-f938a3b204cf
# ╟─e7958ec8-b910-4e03-90bf-f5e50d5af4e7
# ╟─f91c4f28-38b1-4267-8c97-d494f4535cb9
# ╟─77118cea-fd54-4954-944b-c13ed6e0ad1d
# ╠═8be590dc-9b54-4454-88fc-1fa82e91cd5d
# ╠═3427fa7c-41f6-4b08-af59-3eb9d72a7f1a
