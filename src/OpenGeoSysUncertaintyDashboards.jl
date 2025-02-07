module OpenGeoSysUncertaintyDashboards

using OpenGeoSysUncertaintyQuantification
import OpenGeoSysUncertaintyQuantification: StochasticOGS6Parameter
import Makie
import Makie: GridLayout, GridPosition, @L_str, Figure, set_close_to!
import Ogs6InputFileHandler: format_ogs_path
using Printf
using Observables

function standard_dashboard_settings(ogsuqparams::OGSUQParams)
	stoparams = stoch_parameters(ogsuqparams)
	unit_convs = Function[x->x for i = 1:length(stoparams)]
	units = [L"\text{unit}" for i = 1:length(stoparams)]
	mappingfuncs = Any[x->x for i = 1:length(stoparams)]
	settings = Dict{Symbol,Any}(
		:pagewidth => 1200,
		:pageheight => 720,
		:fontsize_main_figure => 20,
		:sliderlinewidth => 10,
		:sliderfontsize => 12f0,
		:draw_freq => 0.5,
		:unit_convs => unit_convs,
		:units => units,
		:mappingfuncs => mappingfuncs,
		:title => "OpenGeoSys UQ Dashboard",
		:plot_frequency => 0.01
	)
end

mutable struct OGSUQDashboard
	figure::Figure
	plot_observables::Dict{Symbol, Observable}
	plot_observable_values::Vector{Any}
	do_plot::Observable{UInt64}
	redraw_timer::Union{Nothing, Timer}
	settings::Dict{Symbol,Any}
	started::Bool
end
sliderobservables(dashboard::OGSUQDashboard) = map(x->dashboard.plot_observables[x], filter(x->occursin("slider_", string(x)), sort(collect(keys(dashboard.plot_observables)))))
sliderobservable_values(dashboard::OGSUQDashboard) = map(x->x[], sliderobservables(dashboard))

function OGSUQDashboard(settings::Dict{Symbol,Any})
	fontsize_mainfig = settings[:fontsize_main_figure]
	title = settings[:title]
	f = Figure(size=(settings[:pagewidth], settings[:pageheight]), fontsize=fontsize_mainfig, title=title)
	plot_observables = Dict{Symbol,Observable}()	
	plot_observable_values = Vector{Any}()	
	return OGSUQDashboard(f, plot_observables, plot_observable_values, Observable{UInt64}(UInt64(0)), nothing, settings, false)
end

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

function parameter_sliders!(dashboard::OGSUQDashboard, pos::GridPosition, stoparams::Vector{StochasticOGS6Parameter}; slider_start_values=zeros(Float64, length(stoparams)))
	settings = dashboard.settings
	sliderlinewidth = settings[:sliderlinewidth]
	mappingfuncs = settings[:mappingfuncs]
	units = settings[:units]
	unit_convs = settings[:unit_convs]
	sliderfontsize = settings[:sliderfontsize]
	sg = Makie.SliderGrid(
		pos,
		[ 
		(
			label = format_ogs_path(stoparam.path), 
			range = -1.0:0.01:1.0, 
			format = value->L"%$(mapintervaltostring(value, mapfun, stoparam, unit_conv)) %$unit", 
			startvalue = 0.0,
			linewidth = sliderlinewidth
			) 
		for (stoparam, mapfun, unit, unit_conv) in zip(stoparams, mappingfuncs, units, unit_convs) 
			]...,
		tellheight=false
	)

	foreach(x->x.fontsize[]=sliderfontsize , sg.labels)
	foreach(x->x.fontsize[]=sliderfontsize , sg.valuelabels)
	foreach(x->x.halign[]=:left , sg.valuelabels)
	foreach((slider,val)->set_close_to!(slider, val), sg.sliders, slider_start_values)
	foreach(tup->dashboard.plot_observables[Symbol("slider_$(tup[1])")] = tup[2].value, enumerate(sg.sliders))
	return nothing
end


function checkplot(dashboard::OGSUQDashboard)
	new_slidervals = map(x->x[], values(dashboard.plot_observables))
	if any(map((x,y)-> x != y, new_slidervals, dashboard.plot_observable_values))
		dashboard.plot_observable_values .= new_slidervals
		dashboard.do_plot[] += 1
	end
	return nothing
end

function OpenGeoSysUncertaintyQuantification.start!(dashboard::OGSUQDashboard)
	dashboard.plot_observable_values = map(x->x[], values(dashboard.plot_observables))
	dashboard.redraw_timer = Timer(cb -> checkplot(dashboard), dashboard.settings[:plot_frequency]; interval=dashboard.settings[:plot_frequency])
	dashboard.started = true
	return dashboard.figure
end

function stop!(dashboard::OGSUQDashboard)
	close(dashboard.redraw_timer)
	sleep(2*dashboard.settings[:plot_frequency])
	dashboard.started = false
end

export standard_dashboard_settings, OGSUQDashboard, parameter_sliders!, start!, stop!, sliderobservables, sliderobservable_values

end # module OpenGeoSysUncertaintyDashboard
