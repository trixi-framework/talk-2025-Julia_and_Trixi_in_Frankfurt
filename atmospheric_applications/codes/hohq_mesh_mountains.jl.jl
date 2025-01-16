using HOHQMesh, GLMakie

MountainProfile = newProject("Mountain_Profile", "Mesh")

setPolynomialOrder!(MountainProfile, 4)
setMeshFileFormat!(MountainProfile, "ABAQUS")

# A background grid is required for the mesh generation.
addBackgroundGrid!(MountainProfile, [250.0*5, 210.0*6, 0.0])

L = 25000.0
H = 21000.0

hc = 2500.0*3
lambda = 4500.0
ac = 10000.0

function funScharMountain(x)

    h = hc*exp(-(x/ac)^2)cos(pi*x/lambda)^2 * (( x + 1000)/7800)^4
    
    t = (x+25000)/50000
    
    return [t x h 0.0]
end

function dome(x)
    t = (x + 25000)/50000
    return [t x sqrt(L^2-x^2) 0.0]
end

x = -25000.0:250:25000
spline_data = zeros(round(length(x)), 4)
spline_data_dome = zeros(round(length(x)), 4)

for i in 1:round(length(x))
    vec = funScharMountain(x[i])
    spline_data[i,:] = vec
    vec = dome(x[length(x) - i + 1])
    spline_data_dome[length(x) - i + 1,:] = vec
end

spline_data[1,3] = 0.0
spline_data[end,3] = 0.0
spline_data_dome[1,3] = 0.0
spline_data_dome[end, 3] = 0.0
spline_data_dome[:,2] = reverse(spline_data_dome[:,2])
outer_spline = newSplineCurve("Bottom", round(length(x)), spline_data)
dome_spline = newSplineCurve("Dome", round(length(x)), spline_data_dome)
addCurveToOuterBoundary!(MountainProfile, outer_spline)
addCurveToOuterBoundary!(MountainProfile, dome_spline)
plotProject!(MountainProfile, MODEL+GRID)
@info "Press enter to generate the mesh and update the plot."
readline()

# Generate the mesh. Saves the mesh file to the directory "out".
generate_mesh(MountainProfile)