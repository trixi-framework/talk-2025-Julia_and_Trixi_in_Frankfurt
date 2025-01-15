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

# ╔═╡ 68362c50-4914-4859-b9c7-1e336bec08e1
begin
	using Pluto
    using VideoIO
	using Images
# Percorso del video
path_to_video = "/home/marco/talk_frankfurt/talk-2025-Julia_and_Trixi_in_Frankfurt/atmospheric_applications/codes/3_rising_bubble/animazione_50_bubble_grid.mp4"

# Apre il video
video = VideoIO.openvideo(path_to_video)

# Funzione per visualizzare il video frame per frame
function show_video(video)
    while !eof(video)  # Continua finché non arriva alla fine del video
        frame = read(video)  # Legge il frame successivo
        display(Images.display(frame))  # Mostra il frame nel notebook forzando il rendering
        sleep(0.1)  # Attende un po' prima di caricare il prossimo frame
    end
end


end

# ╔═╡ ad6ddb85-006a-4ee6-92da-4f56ea4780b0
md""" 
# Atmospheric Applications in Trixi.jl

Some of famous benchmark test cases will be illustrated, showing different possibilities and capabalities of Trixi.jl

- Setting up your own problem
- 
- Different type of Mesh
- How to create Horography
- and more ...
"""

# ╔═╡ 065746f1-1deb-49c4-854d-457a8f7919c2
md""" 
# Setting up your own problem in Trixi.jl

Simulations in Trixi.jl consists of different building blocks
- Set of equations
- Source terms
- Initial conditions
- Boundary Conditions
- Domain size and Mesh (P4est, T8code, StructuredMesh ...)
- Semidiscretization (DGSEM, FDSBP, ...)
- Time Integrator

"""

# ╔═╡ 2f015b91-8c8a-42c7-8264-d4c722591b64
md""" 
# Rising Bubble test case

First we define the equations and the gravity source terms for the momentum and the energy equation.
"""

# ╔═╡ 482f659c-ffbe-44db-bf95-a7bdd86c0640
gamma = 1004.0/717.0

# ╔═╡ e1947159-3084-472c-b9a5-fbdf648042ef
equations  = equations = CompressibleEulerEquations2D(gamma)

# ╔═╡ 012b2045-4172-41f6-98bb-c773b92cf1c3
function source_terms_gravity(u, x, t, equations::CompressibleEulerEquations2D)
    g = 9.81
	rho, _, rho_v2, _ = u
    return SVector(0, 0, -g * rho, -g * rho_v2)
end

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
    potential_temperature = potential_temperature_ref + potential_temperature_perturbation

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

# ╔═╡ 34f5fbd0-abe8-468d-9a62-30d5e5924b97
boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

# ╔═╡ f583fefc-d0a6-4586-8c48-616942d29e98
md""" 
## Defining a simple Mesh

The mesh is a rectangular domain of size $[0, 20 \text{ km}] \text{ x } [0, 10 \text{ km}]$.

The domain is divided into 64 cells along the horizontal direction $x$ and 32 in the vertical direction $z$.
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
                                    boundary_conditions = boundary_conditions)

# ╔═╡ 22a5f545-e62e-4594-92ba-7d7c90d172d4
md""" 
## Creating ODE function and defining the time integrator.

The mesh is a rectangular domain of size $[0, 20 \text{ km}] \text{ x } [0, 10 \text{ km}]$.

The domain is divided into 64 cells along the horizontal direction $x$ and 32 in the vertical direction $z$.
"""

# ╔═╡ f5691db4-fbba-4356-8898-63a33fe56f10
tspan = (0.0, 1000.0)  # 1000 seconds final time

# ╔═╡ bbdb5d15-b207-4b62-a2d6-d496ec8af214
ode = semidiscretize(semi, tspan); nothing

# ╔═╡ d2db4dde-821c-4a4c-80d3-a172b777782a
sol = solve(ode, SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false);

# ╔═╡ f4a2de9f-357c-424f-9619-500d3e524c3f
begin
pd = PlotData2D(sol)
plot(getmesh(pd))
end

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

# ╔═╡ 6a6fd640-396e-41fa-8c6b-977144181ad1
md""" 
## Warped-Curvilinear Mesh.

The mesh is a rectangular domain of size $[0, 20 \text{ km}] \text{ x } [0, 10 \text{ km}]$.

The domain is divided into 64 cells along the horizontal direction $x$ and 32 in the vertical direction $z$.
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

The mesh is a rectangular domain of size $[0, 20 \text{ km}] \text{ x } [0, 10 \text{ km}]$.

The domain is divided into 64 cells along the horizontal direction $x$ and 32 in the vertical direction $z$.
"""

# ╔═╡ 98a46c89-8805-4461-8528-3bd5f4433e31
begin
	trees_per_dimension = (4, 4)
	mesh_amr = P4estMesh(trees_per_dimension,
		polydeg = 3, initial_refinement_level = 2,
		coordinates_min = coordinates_min, coordinates_max = coordinates_max,
		periodicity = (true, false))
end

# ╔═╡ 9051a970-0f5f-4019-b8ee-37b44753865b
begin
	boundary_conditions_amr = Dict(
	:y_neg => boundary_condition_slip_wall,
	:y_pos => boundary_condition_slip_wall)

semi_amr = SemidiscretizationHyperbolic(mesh_amr, 
									equations, 
									initial_condition_rising_bubble, 
									solver,
                                    source_terms = source_terms_gravity,
                                    boundary_conditions = boundary_conditions_amr)
ode_amr = semidiscretize(semi_amr, (0.0, 1000.0))	
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
end

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

# ╔═╡ f3e2f824-dbdc-4748-88a3-ba8ba7e7ddb0
begin
	@inline function cons2perturb(u, equations::CompressibleEulerEquations2D)
	   	c_p = 1004.0
		c_v = 717.0
		p_0 = 100_000.0
	    R = c_p - c_v
		theta_0 = 300.0
	    rho, rho_v1, rho_v2, rho_e = u
	    v1 = rho_v1 / rho
		v2 = rho_v2 / rho
	    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
		rho_theta =  (p / p_0)^(c_v / c_p) * p_0 / R
	    theta = rho_theta/rho - theta_0
	
	    return  rho, v1, v2, theta
	
	end
	
	@inline Trixi.varnames(::typeof(cons2perturb), ::CompressibleEulerEquations2D) = ("rho", "v1", "v2", "theta")
end

# ╔═╡ 92d87248-ba24-45da-9986-0741ffe17082
visualization = VisualizationCallback(interval = 10000,
                                      solution_variables = cons2perturb)

# ╔═╡ f2d35fb1-cba0-419b-b3fc-3fee9a152532
callbacks_amr = CallbackSet(amr_callback)

# ╔═╡ 7ebed9ca-8f57-4971-a274-bb816437642f
begin
sol_amr = solve(ode_amr, SSPRK43(thread = OrdinaryDiffEq.True()),
            maxiters = 1.0e7,
            dt = 1.0,
            save_everystep = false, callback = callbacks_amr);
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

# ╔═╡ 90e6c2a8-9f93-42aa-af6b-1b1c968b98aa
# show video

# ╔═╡ 2e524560-06c7-4dfa-ae8d-d37bcff3214c
show_video(video)

# ╔═╡ b5ac80f0-2fc2-42ea-ad6d-abe43e78cd83

# Chiudi il video
close(video)

# ╔═╡ ef5bba7b-2406-4e82-b509-8cace58bf508
md""" 
## Orography

The mesh is a rectangular domain of size $[0, 20 \text{ km}] \text{ x } [0, 10 \text{ km}]$.

The domain is divided into 64 cells along the horizontal direction $x$ and 32 in the vertical direction $z$.
"""

# ╔═╡ d74493e3-747a-42ac-86ec-3f2b0d01eabb
begin
a = 10000.0
L = 240000.0
H = 30000.0
peak = 2000.0
y_b = peak / (1 + (L/2 / a)^2)
alfa = (H - y_b) * 0.5

f1(s) = SVector(-L/2, y_b + alfa * (s + 1))
f2(s) = SVector(L/2, y_b + alfa * (s + 1))
f3(s) = SVector((s + 1-1) * L/2, peak / (1 + ((s + 1-1) * L/2)^2 / a^2))
f4(s) = SVector((s + 1-1) * L/2, H)

mesh_orography = StructuredMesh((32, 32), (f1, f2, f3, f4), periodicity = (true, false)); nothing
end

# ╔═╡ 5b9eaae3-9278-485d-a119-40004f7bfaec


# ╔═╡ Cell order:
# ╟─ad6ddb85-006a-4ee6-92da-4f56ea4780b0
# ╟─fe03cfcd-bd4b-4c77-99b4-719f51dbe893
# ╟─065746f1-1deb-49c4-854d-457a8f7919c2
# ╟─2f015b91-8c8a-42c7-8264-d4c722591b64
# ╟─482f659c-ffbe-44db-bf95-a7bdd86c0640
# ╠═e1947159-3084-472c-b9a5-fbdf648042ef
# ╠═012b2045-4172-41f6-98bb-c773b92cf1c3
# ╠═61de0650-e3cb-4dfe-823b-65ca76b843dd
# ╠═34f5fbd0-abe8-468d-9a62-30d5e5924b97
# ╟─f583fefc-d0a6-4586-8c48-616942d29e98
# ╟─022ee8fb-c758-4f87-8b60-8103cd8719f5
# ╟─b1290b6e-8a81-41c9-9249-b00d5bfbfd5d
# ╠═4e055027-c7b4-4086-8a8c-a19bf97769c5
# ╠═0e50dafd-a179-4c22-8afa-fa96cd45ff43
# ╠═ff35fbf5-53fd-4e9f-9090-5aab89812b31
# ╠═00f13c5c-9c81-4a68-8885-473d3d61bbf4
# ╠═22a5f545-e62e-4594-92ba-7d7c90d172d4
# ╠═f5691db4-fbba-4356-8898-63a33fe56f10
# ╠═bbdb5d15-b207-4b62-a2d6-d496ec8af214
# ╠═d2db4dde-821c-4a4c-80d3-a172b777782a
# ╟─f4a2de9f-357c-424f-9619-500d3e524c3f
# ╟─54b02548-63e0-4b07-b982-56c51e6f450a
# ╠═e6cdc448-565e-4e50-a726-75847ef13bc5
# ╠═6a6fd640-396e-41fa-8c6b-977144181ad1
# ╠═8fc7fe79-5b78-49a9-9bca-2179460324cc
# ╠═5bb04fa4-36f4-4207-82f1-33a64797ceae
# ╠═956d6c3d-0138-4791-95a9-9f8525570a32
# ╟─0fc2fdc1-4d45-428f-9a47-590409646921
# ╠═98a46c89-8805-4461-8528-3bd5f4433e31
# ╠═9051a970-0f5f-4019-b8ee-37b44753865b
# ╠═28175473-1386-483f-9443-817cad9fd53e
# ╠═498d3262-df89-408e-8ede-f32268894709
# ╠═f3e2f824-dbdc-4748-88a3-ba8ba7e7ddb0
# ╠═92d87248-ba24-45da-9986-0741ffe17082
# ╠═f2d35fb1-cba0-419b-b3fc-3fee9a152532
# ╟─7ebed9ca-8f57-4971-a274-bb816437642f
# ╠═90e6c2a8-9f93-42aa-af6b-1b1c968b98aa
# ╠═68362c50-4914-4859-b9c7-1e336bec08e1
# ╠═2e524560-06c7-4dfa-ae8d-d37bcff3214c
# ╠═b5ac80f0-2fc2-42ea-ad6d-abe43e78cd83
# ╠═ef5bba7b-2406-4e82-b509-8cace58bf508
# ╟─d74493e3-747a-42ac-86ec-3f2b0d01eabb
# ╠═5b9eaae3-9278-485d-a119-40004f7bfaec
