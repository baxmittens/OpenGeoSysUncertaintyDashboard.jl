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

# make a triangle mesh from a quad mesh for visualization only
function make_tri_mesh(faces_quads)
	faces_tri = Matrix{Int64}(undef, 2*size(faces_quads,1), 3)
	for (i,quadface) in enumerate(eachrow(faces_quads))
		j = 2*i-1
		faces_tri[j,:] = quadface[[1,2,3]]
		faces_tri[j+1,:] = quadface[[3,4,1]]
	end
	return faces_tri
end

function mapintervaltostring(x, funn, stoparam, unit_conv)
	mappedval = unit_conv(funn(CPtoStoch(x, stoparam)))
	if abs(mappedval) < 0.001 || abs(mappedval) > 1e4
		return @sprintf("%.3e", mappedval)
	else
		return @sprintf("%.3f", mappedval)
	end
end

xdmf_root_point = DistributedSparseGrids.scaling_weight(first(asg));
asg_ret = deepcopy(xdmf_root_point);
asg_tmp = deepcopy(xdmf_root_point);

vertices = transpose(xdmf_root_point.udata["geometry"][1:2,:,1])
faces_quads = convert(Matrix{Int64},transpose(xdmf_root_point.udata["topology"][1:4,:].+1))
faces = make_tri_mesh(faces_quads)
unit_convs = Function[x->x for i = 1:4]
units = [L"\text{Pa}", L"\text{Pa}", L"\text{Pa}", L"\text{Pa}"]
mappingfuncs = Any[x->x for i = 1:4]

polar = PolarFromCartesian()
r_coords = map(x->polar(x).r, eachrow(vertices))
θ_coords = map(x->polar(x).θ, eachrow(vertices))

const pagewidth = 1200
const pageheight = 720
const fontsize_mainfig = 20
const sliderlinewidth=15
const sliderfontsize=18f0
const draw_freq = 0.5

if @isdefined(redraw_limit) && isopen(redraw_limit)
	close(redraw_limit)
	sleep(2*draw_freq)
end

do_plot = Observable(0)
mainfig=Figure(size=(pagewidth, pageheight), fontsize=fontsize_mainfig, title="Disc with hole");

ax_mesh_displacement_r = Axis(mainfig[1,1], title=L"u_r")
ax_mesh_displacement_θ = Axis(mainfig[1,2], title=L"|u_{\theta}|")
ax_mesh_sigma_rr = Axis(mainfig[3,1], title=L"\sigma_{rr}")
ax_mesh_sigma_θθ = Axis(mainfig[3,2], title=L"\sigma_{\theta\theta}")
ax_mesh_sigma_rθ = Axis(mainfig[3,3], title=L"\sigma_{r\theta}")
ax_mesh_sigma_zz = Axis(mainfig[3,4], title=L"\sigma_{zz}")

allaxes = (ax_mesh_displacement_r, ax_mesh_displacement_θ, ax_mesh_sigma_rr, ax_mesh_sigma_θθ, ax_mesh_sigma_rθ, ax_mesh_sigma_zz)

sg = Makie.SliderGrid(
	mainfig[1,3:4],
	[ 
	(
		label = Ogs6InputFileHandler.format_ogs_path(stoparam.path), 
		range = -1.0:0.01:1.0, 
		format = value->L"%$(mapintervaltostring(value, mapfun, stoparam, unit_conv)) %$unit", 
		startvalue = 0.0,
		linewidth = sliderlinewidth
		) 
	for (stoparam, mapfun, unit, unit_conv) in zip(stoch_parameters(ogsuqasg), mappingfuncs, units, unit_convs) 
		]...,
	tellheight=false
)

sg_options = Makie.SliderGrid(
	mainfig[2,3:4],
	(
		label = "displ. mult.", 
		range = 1:1:10, 
		startvalue = 1,
		linewidth = sliderlinewidth
	),
	tellheight=false
)

foreach(x->x.fontsize[]=sliderfontsize , sg.labels)
foreach(x->x.fontsize[]=sliderfontsize , sg.valuelabels)
foreach(x->x.halign[]=:left , sg.valuelabels)
sliderobservables = [s.value for s in sg.sliders]
set_close_to!(sg.sliders[end], 1.0)

old_slidervals = ones(Float64, length(sliderobservables))
function checkplot(sliderobservables=sliderobservables, old_slidervals=old_slidervals)
	new_slidervals = map(x->x[], sliderobservables)
	if any(map((x,y)->!isapprox(x,y,atol=0.01), new_slidervals, old_slidervals))
		old_slidervals .= new_slidervals
		do_plot[] += 1
	end
	return nothing
end

xdmf_asg = map!(Observable{Any}(), do_plot) do _plotnow_
	x = map(x->x[], sliderobservables)
	interpolate!(asg_ret, asg_tmp, asg, x)
	return asg_ret
end

displacement_r = map!(Observable{Any}(), xdmf_asg) do xdmf
	return map(x->polar(x).r, eachcol(xdmf["displacement"][1:2,:,end]))
end
displacement_r_limits = map!(Observable{Any}(), displacement_r) do u
	lims = minimum(u),maximum(u)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end
displacement_θ = map!(Observable{Any}(), xdmf_asg) do xdmf
	return map(x->abs(polar(x).θ), eachcol(xdmf["displacement"][1:2,:,end]))
end
displacement_θ_limits = map!(Observable{Any}(), displacement_θ) do u
	lims = minimum(u),maximum(u)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end
sigma_rr = map!(Observable{Any}(), xdmf_asg) do xdmf
	return xdmf["sigma"][1,:,end] .* cos.(θ_coords).^2.0 .+ xdmf["sigma"][2,:,end] .* sin.(θ_coords).^2.0 .+ xdmf["sigma"][4,:,end] .* sin.(2.0.*θ_coords)
end
sigma_rr_limits = map!(Observable{Any}(), sigma_rr) do σ
	lims = minimum(σ),maximum(σ)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end
sigma_θθ = map!(Observable{Any}(), xdmf_asg) do xdmf
	return xdmf["sigma"][1,:,end] .* sin.(θ_coords).^2.0 .+ xdmf["sigma"][2,:,end] .* cos.(θ_coords).^2.0 .- xdmf["sigma"][4,:,end] .* sin.(2.0*θ_coords)
end
sigma_θθ_limits = map!(Observable{Any}(), sigma_θθ) do σ
	lims = minimum(σ),maximum(σ)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end
sigma_rθ = map!(Observable{Any}(), xdmf_asg) do xdmf
	return sin.(θ_coords) .* cos.(θ_coords) .* (xdmf["sigma"][2,:,end] .- xdmf["sigma"][1,:,end]) .+ xdmf["sigma"][4,:,end] .* cos.(2.0*θ_coords)
end
sigma_rθ_limits = map!(Observable{Any}(), sigma_rθ) do σ
	lims = minimum(σ),maximum(σ)
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end
sigma_zz = map!(Observable{Any}(), xdmf_asg) do xdmf
	return xdmf["sigma"][3,:,end]
end
sigma_zz_limits = map!(Observable{Any}(), xdmf_asg) do xdmf
	lims = minimum(xdmf["sigma"][3,:,end]),maximum(xdmf["sigma"][3,:,end])
	if abs(lims[2]-lims[1]) < 1e-6
		return -1.0,1.0
	else
		return lims
	end
end

plot_vertices_x = map!(Observable{Any}(), xdmf_asg, sg_options.sliders[1].value) do xdmf,mult
	verts = vertices[:,1] .+ xdmf["displacement"][1,:,end].*mult
	foreach(x->xlims!(x, (minimum(verts), maximum(verts))), allaxes)	
	return verts
end
plot_vertices_y = map!(Observable{Any}(), xdmf_asg, sg_options.sliders[1].value) do xdmf,mult
	verts = vertices[:,2] .+ xdmf["displacement"][2,:,end].*mult
	foreach(x->ylims!(x, (minimum(verts), maximum(verts))), allaxes)	
	return verts
end

tricontourf!(ax_mesh_displacement_r, plot_vertices_x, plot_vertices_y, displacement_r, triangulation = faces)
tricontourf!(ax_mesh_displacement_θ, plot_vertices_x, plot_vertices_y, displacement_θ, triangulation = faces)
tricontourf!(ax_mesh_sigma_rr, plot_vertices_x, plot_vertices_y, sigma_rr, triangulation = faces)
tricontourf!(ax_mesh_sigma_θθ, plot_vertices_x, plot_vertices_y, sigma_θθ, triangulation = faces)
tricontourf!(ax_mesh_sigma_rθ, plot_vertices_x, plot_vertices_y, sigma_rθ, triangulation = faces)
tricontourf!(ax_mesh_sigma_zz, plot_vertices_x, plot_vertices_y, sigma_zz, triangulation = faces)

tickformat=values->[@sprintf("%.3f", val) for val in values]
Colorbar(mainfig[2,1], limits=displacement_r_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(mainfig[2,2], limits=displacement_θ_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(mainfig[4,1], limits=sigma_rr_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(mainfig[4,2], limits=sigma_θθ_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(mainfig[4,3], limits=sigma_rθ_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)
Colorbar(mainfig[4,4], limits=sigma_zz_limits, vertical=false, flipaxis=false, tickformat=tickformat, ticklabelsize=10)

redraw_limit = Timer(cb -> checkplot(), 0.1; interval=0.5)

display(mainfig)