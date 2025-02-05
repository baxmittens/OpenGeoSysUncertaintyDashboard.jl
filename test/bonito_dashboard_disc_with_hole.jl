using OpenGeoSysUncertaintyQuantification
using DistributedSparseGrids
using Ogs6InputFileHandler
using WGLMakie
using Printf
using CoordinateTransformations
using Bonito

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
units = [L"\text{N}", L"\text{N}", L"\text{N}", L"\text{N}"]
mappingfuncs = Any[x->x for i = 1:4]

polar = PolarFromCartesian()
r_coords = map(x->polar(x).r, eachrow(vertices))
θ_coords = map(x->polar(x).θ, eachrow(vertices))

const pagewidth = 1200
const pageheight = 720
const fontsize_mainfig = 20
const sliderlinewidth=10
const sliderfontsize=12f0
const draw_freq = 0.5

if @isdefined(redraw_limit) && isopen(redraw_limit)
	close(redraw_limit)
	sleep(2*draw_freq)
end

check_plot = Observable(0)

app = App() do session::Session

	do_plot = Observable(0)
	mainfig=Figure(size=(pagewidth, pageheight), fontsize=fontsize_mainfig, title="Disc with hole");
	#gridcontrols = mainfig[1:2,3:4] = GridLayout()
	
	#rowgap!(gridcontrols,10)
	
	ax_bc = Axis(mainfig[1,3:4], aspect = DataAspect())
	ax_mesh_displacement_r = Axis(mainfig[1,1], title=L"u_r")
	ax_mesh_displacement_θ = Axis(mainfig[1,2], title=L"|u_{\theta}|")
	ax_mesh_sigma_rr = Axis(mainfig[3,1], title=L"\sigma_{rr}")
	ax_mesh_sigma_θθ = Axis(mainfig[3,2], title=L"\sigma_{\theta\theta}")
	ax_mesh_sigma_rθ = Axis(mainfig[3,3], title=L"\sigma_{r\theta}")
	ax_mesh_sigma_zz = Axis(mainfig[3,4], title=L"\sigma_{zz}")
	
	allaxes = (ax_mesh_displacement_r, ax_mesh_displacement_θ, ax_mesh_sigma_rr, ax_mesh_sigma_θθ, ax_mesh_sigma_rθ, ax_mesh_sigma_zz)
	
	hidedecorations!(ax_bc)
	hidespines!(ax_bc)
	
	#sg = Makie.SliderGrid(
	#	gridcontrols[4,1],
	#	[ 
	#	(
	#		label = Ogs6InputFileHandler.format_ogs_path(stoparam.path), 
	#		range = -1.0:0.01:1.0, 
	#		format = value->L"%$(mapintervaltostring(value, mapfun, stoparam, unit_conv)) %$unit", 
	#		startvalue = 0.0,
	#		linewidth = sliderlinewidth
	#		) 
	#	for (stoparam, mapfun, unit, unit_conv) in zip(stoch_parameters(ogsuqasg), mappingfuncs, units, unit_convs) 
	#		]...,
	#	tellheight=false
	#)
	#
	#sg_options = Makie.SliderGrid(
	#	gridcontrols[6,1],
	#	(
	#		label = "displ. mult.", 
	#		range = 1:1:10, 
	#		startvalue = 1,
	#		linewidth = sliderlinewidth
	#	),
	#	tellheight=false
	#)

	bonito_sliders = [i == 4 ? Bonito.Slider(-1.0:0.01:1.0, value=1.0) : Bonito.Slider(-1.0:0.01:1.0, value=0.0) for i = 1:4]
	
	bonito_control_slider = Bonito.Slider(1:1:10, value=1)
	#foreach(x->x.fontsize[]=sliderfontsize , sg.labels)
	#foreach(x->x.fontsize[]=sliderfontsize , sg.valuelabels)
	#foreach(x->x.halign[]=:left , sg.valuelabels)
	#foreach(x->x.fontsize[]=sliderfontsize , sg_options.labels)
	#foreach(x->x.fontsize[]=sliderfontsize , sg_options.valuelabels)
	#foreach(x->x.halign[]=:left , sg_options.valuelabels)
	
	sliderobservables = [s.value for s in bonito_sliders]
	#set_close_to!(sg.sliders[end], 1.0)
	
	old_slidervals = ones(Float64, length(sliderobservables))
	old_multval = [0]
	_ = map!(Observable{Any}(), check_plot) do _checknow_
		new_slidervals = map(x->x[], sliderobservables)
		new_multval = bonito_control_slider.value[]
		if any(map((x,y)->!isapprox(x,y,atol=0.01), new_slidervals, old_slidervals)) || old_multval[1] != new_multval
			old_slidervals .= new_slidervals
			old_multval[1] = new_multval
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
	
	plot_vertices_x = map!(Observable{Any}(), xdmf_asg) do xdmf
		mult = bonito_control_slider.value[]
		verts = vertices[:,1] .+ xdmf["displacement"][1,:,end].*mult
		foreach(x->xlims!(x, (minimum(verts), maximum(verts))), allaxes)	
		return verts
	end
	plot_vertices_y = map!(Observable{Any}(), xdmf_asg) do xdmf
		mult = bonito_control_slider.value[]
		verts = vertices[:,2] .+ xdmf["displacement"][2,:,end].*mult
		foreach(x->ylims!(x, (minimum(verts), maximum(verts))), allaxes)	
		return verts
	end
	
	rightarrowx = map!(Observable{Any}(), do_plot) do _plotnow_
		scalex = sliderobservables[1][]
		_rightarrowx = [scalex > 0.0 ? 10.0 : 10 - 2*scalex for i = 1:5]
		return _rightarrowx
	end
	rightarrowdx = map!(Observable{Any}(), do_plot) do _plotnow_
		scalex = sliderobservables[1][]
		_rightarrowdx = [2*scalex for i = 1:5]
		return _rightarrowdx
	end
	rightarrowx_visible = map!(Observable{Any}(), do_plot) do _plotnow_
		scalex = sliderobservables[1][]
		if isapprox(scalex, 0.0, atol=0.1)
			return false
		end
		return true
	end
	rightarrowdy = map!(Observable{Any}(), do_plot) do _plotnow_
		scaley = sliderobservables[2][]
		_rightarrowdy = [2*scaley for i = 1:5]
		return _rightarrowdy
	end
	rightarrowy_visible = map!(Observable{Any}(), do_plot) do _plotnow_
		scaley = sliderobservables[2][]
		if isapprox(scaley, 0.0, atol=0.1)
			return false
		end
		return true
	end
	
	toparrowdx = map!(Observable{Any}(), do_plot) do _plotnow_
		scalex = sliderobservables[3][]
		_toparrowdx = [2*scalex for i = 1:5]
		return _toparrowdx
	end
	toparrowx_visible = map!(Observable{Any}(), do_plot) do _plotnow_
		scalex = sliderobservables[3][]
		if isapprox(scalex, 0.0, atol=0.1)
			return false
		end
		return true
	end
	
	toparrowy = map!(Observable{Any}(), do_plot) do _plotnow_
		scaley = sliderobservables[4][]
		_toparrowy = [scaley > 0.0 ? 10.0 : 10 - 2*scaley for i = 1:5]
		return _toparrowy
	end
	toparrowdy = map!(Observable{Any}(), do_plot) do _plotnow_
		scaley = sliderobservables[4][]
		_toparrowdy = [2*scaley for i = 1:5]
		return _toparrowdy
	end
	toparrowy_visible = map!(Observable{Any}(), do_plot) do _plotnow_
		scalex = sliderobservables[4][]
		if isapprox(scalex, 0.0, atol=0.1)
			return false
		end
		return true
	end
	
	ylims!(ax_bc, (-2, 13))
	xlims!(ax_bc, (-2, 13))
	
	mesh!(ax_bc, vertices, faces)
	
	arrows!(ax_bc, rightarrowx, [i for i = 1.0:2.0:9.0], rightarrowdx, [0.0 for i = 1:5], visible=rightarrowx_visible)
	arrows!(ax_bc, [10.5 for i = 1:5], [i for i = 1.0:2.0:9.0], [0.0 for i = 1:5], rightarrowdy, color=:blue, visible=rightarrowy_visible)
	arrows!(ax_bc, [i for i = 1.0:2.0:9.0], [10.5 for i = 1:5], toparrowdx, [0.0 for i = 1:5],  color=:blue, visible=toparrowx_visible)
	arrows!(ax_bc, [i for i = 1.0:2.0:9.0], toparrowy, [0.0 for i = 1:5], toparrowdy, visible=toparrowy_visible)
	
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
	return DOM.div(DOM.div(mainfig), DOM.div(map(x->DOM.div(x, x.value), bonito_sliders), DOM.div(bonito_control_slider,bonito_control_slider.value)))
end

redraw_limit = Timer(cb -> check_plot[]+=1 , 1.0; interval=1.0)

display(app)