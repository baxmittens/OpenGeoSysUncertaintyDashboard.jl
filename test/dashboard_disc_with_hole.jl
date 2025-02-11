using OpenGeoSysUncertaintyDashboards
using OpenGeoSysUncertaintyQuantification
using DistributedSparseGrids
using Ogs6InputFileHandler
using GLMakie
using Printf
using CoordinateTransformations

# generate/load the stochastic model
ogsuqparams = OGSUQParams("../projects/DiscWithHole/StochasticOGSModelParams.xml", "../projects/DiscWithHole/SampleMethodParams.xml")
ogsuqasg = init(ogsuqparams)
start!(ogsuqasg)
asg = ogsuqasg.asg

# instantiate the dashboard
settings = standard_dashboard_settings(ogsuqparams)
dashboard = OGSUQDashboard(settings)

# get a prototype of the result file/topology for efficient sparse grid interpolation
xdmf_root_point = DistributedSparseGrids.scaling_weight(first(asg));
asg_ret = deepcopy(xdmf_root_point);
asg_tmp = deepcopy(xdmf_root_point);

# get vertices and faces
vertices = transpose(xdmf_root_point.udata["geometry"][1:2,:,1])
faces_quads = convert(Matrix{Int64},transpose(xdmf_root_point.udata["topology"][1:4,:].+1))
faces = OpenGeoSysUncertaintyDashboards.make_tri_mesh(faces_quads)

# get Polar coordinates for some additional result postprocessing
polar = PolarFromCartesian()
r_coords = map(x->polar(x).r, eachrow(vertices))
θ_coords = map(x->polar(x).θ, eachrow(vertices))

# set up the main layout of the dashboard
gridcontrols = dashboard.figure[1:2,3:4] = GridLayout()
ax_bc = Axis(gridcontrols[1:2,1], aspect = DataAspect())
ax_mesh_displacement_r = Axis(dashboard.figure[1,1], title=L"u_r")
ax_mesh_displacement_θ = Axis(dashboard.figure[1,2], title=L"|u_{\theta}|")
ax_mesh_sigma_rr = Axis(dashboard.figure[3,1], title=L"\sigma_{rr}")
ax_mesh_sigma_θθ = Axis(dashboard.figure[3,2], title=L"\sigma_{\theta\theta}")
ax_mesh_sigma_rθ = Axis(dashboard.figure[3,3], title=L"\sigma_{r\theta}")
ax_mesh_sigma_zz = Axis(dashboard.figure[3,4], title=L"\sigma_{zz}")
# bundle all axes for rescaling axes if result displacements are amplified
allaxes = (ax_mesh_displacement_r, ax_mesh_displacement_θ, ax_mesh_sigma_rr, ax_mesh_sigma_θθ, ax_mesh_sigma_rθ, ax_mesh_sigma_zz)
# hide all decorations for the plot showing the boundary conditions
hidedecorations!(ax_bc)
hidespines!(ax_bc)

# set the start positions of the sliders
slider_start_values = [0.0, 0.0, 0.0, 0.5]
# initialized the sliders with the dashboard
parameter_sliders!(dashboard, gridcontrols[3,1], stoch_parameters(ogsuqasg), slider_start_values=slider_start_values)

# set up an additional slider for displacement amplification
additional_slider = Makie.SliderGrid(
	gridcontrols[4,1],
	(
		label = "displ. mult.", 
		range = 1:1:10, 
		startvalue = 1,
		linewidth = dashboard.settings[:sliderlinewidth]
	),
	tellheight=false
)
foreach(x->x.fontsize[]=dashboard.settings[:sliderfontsize] , additional_slider.labels)
foreach(x->x.fontsize[]=dashboard.settings[:sliderfontsize] , additional_slider.valuelabels)
foreach(x->x.halign[]=:left , additional_slider.valuelabels)

# since this should update the plot, this has to be added to the plot_obserables
dashboard.plot_observables[:additional_slider] = first(additional_slider.sliders).value

## Creation of Observables for reactive plotting

# here, the interpolation with the sparse grid takes place.
# the result is an XDMF3File
# the frequency with this event can repeat itself is limited by `dashboard.settings[:plot_frequency]
xdmf_asg = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
	x = sliderobservable_values(dashboard)
	interpolate!(asg_ret, asg_tmp, asg, x)
	return asg_ret
end

# displacement field r in polar coordinates for colors in tricontourf plot
displacement_r = map!(Observable{Any}(), xdmf_asg) do xdmf
	return map(x->polar(x).r, eachcol(xdmf["displacement"][1:2,:,end]))
end

# displacement field r limits for colorbar
displacement_r_limits = map!(Observable{Any}(), displacement_r) do u
	lims = minimum(u),maximum(u)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end

# displacement field θ in polar coordinates for colors in tricontourf plot
displacement_θ = map!(Observable{Any}(), xdmf_asg) do xdmf
	return map(x->abs(polar(x).θ), eachcol(xdmf["displacement"][1:2,:,end]))
end

# displacement field θ limits for colorbar
displacement_θ_limits = map!(Observable{Any}(), displacement_θ) do u
	lims = minimum(u),maximum(u)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end

# sigma rr θ in polar coordinates for colors in tricontourf plot
sigma_rr = map!(Observable{Any}(), xdmf_asg) do xdmf
	return xdmf["sigma"][1,:,end] .* cos.(θ_coords).^2.0 .+ xdmf["sigma"][2,:,end] .* sin.(θ_coords).^2.0 .+ xdmf["sigma"][4,:,end] .* sin.(2.0.*θ_coords)
end

# sigma rr limits for colorbar
sigma_rr_limits = map!(Observable{Any}(), sigma_rr) do σ
	lims = minimum(σ),maximum(σ)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end

# sigma θθ in polar coordinates for colors in tricontourf plot
sigma_θθ = map!(Observable{Any}(), xdmf_asg) do xdmf
	return xdmf["sigma"][1,:,end] .* sin.(θ_coords).^2.0 .+ xdmf["sigma"][2,:,end] .* cos.(θ_coords).^2.0 .- xdmf["sigma"][4,:,end] .* sin.(2.0*θ_coords)
end

# sigma θθ limits for colorbar
sigma_θθ_limits = map!(Observable{Any}(), sigma_θθ) do σ
	lims = minimum(σ),maximum(σ)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end

# sigma rθ in polar coordinates for colors in tricontourf plot
sigma_rθ = map!(Observable{Any}(), xdmf_asg) do xdmf
	return sin.(θ_coords) .* cos.(θ_coords) .* (xdmf["sigma"][2,:,end] .- xdmf["sigma"][1,:,end]) .+ xdmf["sigma"][4,:,end] .* cos.(2.0*θ_coords)
end

# sigma rθ limits for colorbar
sigma_rθ_limits = map!(Observable{Any}(), sigma_rθ) do σ
	lims = minimum(σ),maximum(σ)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end

# sigma zz for colors in tricontourf plot
sigma_zz = map!(Observable{Any}(), xdmf_asg) do xdmf
	return xdmf["sigma"][3,:,end]
end

# sigma zz limits for colorbar
sigma_zz_limits = map!(Observable{Any}(), xdmf_asg) do xdmf
	lims = minimum(xdmf["sigma"][3,:,end]),maximum(xdmf["sigma"][3,:,end])
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end

# vertex coordinates x for all plots
# gets amplified with additional slider
plot_vertices_x = map!(Observable{Any}(), xdmf_asg) do xdmf
	mult = dashboard.plot_observables[:additional_slider][]
	verts = vertices[:,1] .+ xdmf["displacement"][1,:,end].*mult
	foreach(x->xlims!(x, (minimum(verts), maximum(verts))), allaxes)	
	return verts
end

# vertex coordinates y for all plots
# gets amplified with additional slider
plot_vertices_y = map!(Observable{Any}(), xdmf_asg) do xdmf
	mult = dashboard.plot_observables[:additional_slider][]
	verts = vertices[:,2] .+ xdmf["displacement"][2,:,end].*mult
	foreach(x->ylims!(x, (minimum(verts), maximum(verts))), allaxes)	
	return verts
end

## the following block are all Observables for reactive representation of boundary conditions
begin
	rightarrowx = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scalex = x[1]
		rightarrowx = [scalex > 0.0 ? 10.0 : 10 - 2*scalex for i = 1:5]
		return rightarrowx
	end
	rightarrowdx = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scalex = x[1]
		rightarrowdx = [2*scalex for i = 1:5]
		return rightarrowdx
	end
	rightarrowx_visible = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scalex = x[1]
		if isapprox(scalex, 0.0, atol=0.05)
			return false
		end
		return true
	end
	rightarrowdy = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scaley = x[2]
		rightarrowdx = [2*scaley for i = 1:5]
		return rightarrowdx
	end
	rightarrowy_visible = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scaley = x[2]
		if isapprox(scaley, 0.0, atol=0.05)
			return false
		end
		return true
	end
	toparrowdx = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scalex = x[3]
		rightarrowdx = [2*scalex for i = 1:5]
		return rightarrowdx
	end
	toparrowx_visible = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scalex = x[3]
		if isapprox(scalex, 0.0, atol=0.05)
			return false
		end
		return true
	end
	toparrowy = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scaley = x[4]
		rightarrowy = [scaley > 0.0 ? 10.0 : 10 - 2*scaley for i = 1:5]
		return rightarrowy
	end
	toparrowdy = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scaley = x[4]
		rightarrowdy = [2*scaley for i = 1:5]
		return rightarrowdy
	end
	toparrowy_visible = map!(Observable{Any}(), dashboard.do_plot) do _plotnow_
		x = sliderobservable_values(dashboard)
		scalex = x[4]
		if isapprox(scalex, 0.0, atol=0.05)
			return false
		end
		return true
	end
end

# set the limits for boundary condition plot
ylims!(ax_bc, (-2, 13))
xlims!(ax_bc, (-2, 13))
# plot the boundary condition plot
mesh!(ax_bc, vertices, faces)
arrows!(ax_bc, rightarrowx, [i for i = 1.0:2.0:9.0], rightarrowdx, [0.0 for i = 1:5], visible=rightarrowx_visible, arrowsize=5)
arrows!(ax_bc, [10.5 for i = 1:5], [i for i = 1.0:2.0:9.0], [0.0 for i = 1:5], rightarrowdy, visible=rightarrowy_visible, arrowsize=5, color=:blue)
arrows!(ax_bc, [i for i = 1.0:2.0:9.0], [10.5 for i = 1:5], toparrowdx, [0.0 for i = 1:5], visible=toparrowx_visible, arrowsize=5, color=:blue)
arrows!(ax_bc, [i for i = 1.0:2.0:9.0], toparrowy, [0.0 for i = 1:5], toparrowdy, visible=toparrowy_visible, arrowsize=5)

# plot all result plots
tricontourf!(ax_mesh_displacement_r, plot_vertices_x, plot_vertices_y, displacement_r, triangulation = faces)
tricontourf!(ax_mesh_displacement_θ, plot_vertices_x, plot_vertices_y, displacement_θ, triangulation = faces)
tricontourf!(ax_mesh_sigma_rr, plot_vertices_x, plot_vertices_y, sigma_rr, triangulation = faces)
tricontourf!(ax_mesh_sigma_θθ, plot_vertices_x, plot_vertices_y, sigma_θθ, triangulation = faces)
tricontourf!(ax_mesh_sigma_rθ, plot_vertices_x, plot_vertices_y, sigma_rθ, triangulation = faces)
tricontourf!(ax_mesh_sigma_zz, plot_vertices_x, plot_vertices_y, sigma_zz, triangulation = faces)

# plot all colorbars
tickformat=values->[@sprintf("%.3f", val) for val in values]
Colorbar(dashboard.figure[2,1], limits=displacement_r_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(dashboard.figure[2,2], limits=displacement_θ_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(dashboard.figure[4,1], limits=sigma_rr_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(dashboard.figure[4,2], limits=sigma_θθ_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(dashboard.figure[4,3], limits=sigma_rθ_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(dashboard.figure[4,4], limits=sigma_zz_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)

start!(dashboard)
