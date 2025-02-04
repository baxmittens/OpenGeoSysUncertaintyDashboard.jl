using OpenGeoSysUncertaintyQuantification

PATH = @__DIR__

projectfile= joinpath(PATH, "RefPrj/disc_with_hole.prj")
pathes = generatePossibleStochasticParameters(projectfile, joinpath(PATH, "StochasticParameters.xml"))
inds = findall(x->occursin("?neumann_force", x), pathes)
writeStochasticParameters(pathes[inds], joinpath(PATH, "StochasticParameters.xml"))

simcall = OpenGeoSysUncertaintyQuantification.install_ogs() #needs python
additionalprojecfilespath=joinpath(PATH, "Misc")
outputpath= joinpath(PATH, "Res")
postprocfiles=["disc_with_hole_disc_with_hole.xdmf"]
stochmethod=AdaptiveHierarchicalSparseGrid
n_workers = 20

stochparampathes = loadStochasticParameters(joinpath(PATH,"StochasticParameters.xml"))
	
stochasticmodelparams = generateStochasticOGSModell(
	projectfile,
	simcall,
	additionalprojecfilespath,
	postprocfiles,
	stochparampathes,
	outputpath,
	stochmethod,
	n_workers,
	joinpath(PATH,"user_functions.jl"),
	joinpath(PATH,"StochasticOGSModelParams.xml")
	)

using Distributions

for stop in stochasticmodelparams.stochparams 
	stop.lower_bound = -0.01
	stop.upper_bound = 0.01
	stop.dist = Uniform(-0.01,0.01)
end
write(stochasticmodelparams)
display(stochasticmodelparams)
samplemethodparams = generateSampleMethodModel(stochasticmodelparams, joinpath(PATH, "SampleMethodParams.xml"))
display(samplemethodparams)
