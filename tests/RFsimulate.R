library("geostatsp")

Sys.setenv(
  OMP_NUM_THREADS = "2",
  OPENBLAS_NUM_THREADS = "2",
  MKL_NUM_THREADS = "2",
  BLIS_NUM_THREADS = "2",
  VECLIB_MAXIMUM_THREADS = "2",   # macOS Accelerate
  NUMEXPR_NUM_THREADS = "2"
)
options(mc.cores = 2)

model <- c(var=5, range=20,shape=0.5)

# any old crs
theCrs = "+proj=utm +zone=17 +datum=NAD27 +units=m +no_defs"

# don't test using the randomFields package, it's currently broken on some R builds
options(useRandomFields = FALSE)

  myraster = rast(nrows=20,ncols=20,extent = ext(100,110,100,110), 
    crs=theCrs)

set.seed(0)
simu = RFsimulate(model = rbind(a=model, b=model+0.1), 
  x=myraster, n=4
)


par(mfrow=c(ceiling(length(names(simu))/2), 2))

for(D in 1:length(names(simu))) {
  plot(simu[[D]])
}


xPoints = suppressWarnings(
  as.points(myraster)[
      sample(ncell(myraster), 12),]
)
  simu2 = RFsimulate(
    model = rbind(a=model, b=model+0.1), 
    x= xPoints
   )
  
  
    par(mfrow=c(nlyr(simu),2))
    for(D in 1:nlyr(simu)) {
      plot(simu[[D]])
      plot(simu2)
    }
  
  
  
  data("swissRain")
  swissRain = unwrap(swissRain)
  swissAltitude = unwrap(swissAltitude)
  swissBorder = unwrap(swissBorder)
  swissRain$sqrtrain = sqrt(swissRain$rain)
  
# estimate parameters
  
  
# isotropic
  swissRes =  lgm(data=swissRain, 
    grid=20, formula="sqrtrain",
    covariates=swissAltitude,   
    shape=1, fixShape=TRUE,
    aniso=FALSE, fixNugget=FALSE,
    nuggetInPrediction=FALSE
  )
  
  
  # anisotropic
  swissRes =  lgm("sqrtrain",
    swissRain, grid=20, 
    covariates=swissAltitude,   
    shape=1, fixShape=TRUE,
    aniso=TRUE, fixNugget=FALSE,
    nuggetInPrediction=FALSE
  )
  
  
  
# uncoinditional simulation
  swissSim = RFsimulate(
    model=swissRes$param,
    x=swissRes$predict,
    n=3
  )
  
  
# simulate from the random effect conditional on
#   the observed data
  
  # there's a bug in RandomFields and this won't run
  options(useRandomFields=FALSE)
  swissSim = RFsimulate(
    model=swissRes$param,
    data=swissRes$data[,'resid'],
    x=swissRes$predict,
    err.model=swissRes$param["nugget"]+1,
    n=3
  )
  
# plot the simulated random effect
  plot(swissSim[[1]])
  plot(swissBorder, add=TRUE)
  
# now with multiple parameter sets 
  swissSim = RFsimulate(model=
      rbind(
        swissRes$param,
        swissRes$param*0.9),
    data=swissRes$data[,'resid'],
    x=swissRes$predict,
    err.model=c(1, 0.9)*swissRes$param["nugget"]+1,
    n=3
  )
# plot the simulated random effect
  plot(swissSim[[1]])
  plot(swissBorder, add=TRUE)
  
# and multiple simulations
# now with multiple datasets 
  swissSim = RFsimulate(model=
      rbind(
        swissRes$param,	
        0.99*swissRes$param,
        1.01*swissRes$param),
    data=swissRes$data[,rep('resid',3)],
    err.model=c(1, 0.99, 1.01)*swissRes$param["nugget"]+1,
    x=swissRes$predict,
    n=3
  )
# plot the simulated random effect
  plot(swissSim[[1]])
  plot(swissBorder, add=TRUE)
  
  
  
# create a small raster of elevation data
  swissAltSmall = aggregate(swissAltitude,fact=5)
  swissAltSmall = resample(swissAltSmall, swissSim)
  
# calculate the fixed effects portion of the rainfall process
  rainMean = swissRes$param["(Intercept)"] +
    swissRes$param[ "CHE_alt" ] * swissAltSmall
  
# define a function to identify the location of maximum rainfall	
  
  swissLocation = terra::global(swissSim,   which.max)
  swissLocation = xyFromCell(swissSim, unlist(swissLocation))
  plot(swissRes$predict[["predict"]])
  plot(swissBorder, add=TRUE)
  points(swissLocation)
  
  


