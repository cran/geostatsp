
library('geostatsp')

grd <- rast(extent = c(0,10,0,10), res = c(1,1))
polys <- as.polygons(grd)
centroids <- crds(as.points(polys))[seq(from=1, by=4, len=ncell(grd)),]
x <- centroids[,1]
y <- centroids[,2]
z <- 1.4 + 0.1*x + 0.2*y + 0.002*x*x
xpoly <- polys
values(xpoly) =data.frame(x=x, y=y, z=z, row.names=row.names(polys))

names(xpoly)=paste("stuff", 
    1:length(names(xpoly)), sep="")


template = squareRaster(xpoly, 100)

thebrick = spdfToBrick(
    x=xpoly,
    template=template,
    pattern='^stuff[[:digit:]]$'
    )
    
plot(thebrick[[1]])   

