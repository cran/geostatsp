### sample usage:
###  lancs =  getTiles(c(-2.842,-2.7579),c(54.0295,54.063),12,path="http://tile.openstreetmap.org/",maxTiles=60,verbose=TRUE)
###

.tile2boundingBox <- function(x,y,zoom){

  n = .tile2lat(y,zoom)
  s = .tile2lat(y+1,zoom)
  w = .tile2lon(x,zoom)
  e = .tile2lon(x+1,zoom)
  return(c(s,w,n,e))
  
}

.tile2extentMercator <- function(I, J, zoom) {
	
	# formula from http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
	n = 2 ^ zoom
	lon_rad = c(I, I+1) / n * 2*pi - pi
	lon_deg = lon_rad *180/pi
	lat_rad = atan(sinh(pi * (1 - 2 * c(J+1,J) / n)))
	lat_deg = lat_rad * 180.0 / pi

	thePoints = SpatialPoints(cbind(lon_deg, lat_deg), CRS("+proj=longlat"))
	thePointsMerc = spTransform(thePoints, 
			CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
			)

	round(extent(thePointsMerc))			
}


.tile2lat <- function(y,zoom){
  n = pi - ((2.0 * pi * y) / (2.0^zoom))
  return ( 180.0 / pi * atan(0.5 * (exp(n) - exp(-n))))
}

.tile2lon <- function(x,zoom){
  return(((x/(2^zoom)*360)-180))
}


.lonlat2tile <- function(lon,lat,zoom){
  xtile = floor(((lon + 180) / 360) * 2^zoom)
  ytile = floor((1 - log(tan(lat*pi/180) + 1 / cos(lat*pi/180)) / pi) /2 * 2^zoom)
  return(c(xtile,ytile))
}

.getTileBounds <- function(xlim,ylim,zoom){
  LL = .lonlat2tile(xlim[1],ylim[1],zoom)
  UR = .lonlat2tile(xlim[2],ylim[2],zoom)
  return(list(LL,UR))
}

nTiles <- function(xlim,ylim,zoom){
  tb = .getTileBounds(xlim,ylim,zoom)
  nt = (tb[[2]][1]-tb[[1]][1]+1)*(tb[[1]][2]-tb[[2]][2]+1)
  return(nt)
}

getTilePaths <- function(xlim,ylim,zoom,path){
  tileBounds = .getTileBounds(xlim,ylim,zoom)
  LL = tileBounds[[1]]
  UR = tileBounds[[2]]
  tileData = list()
  i = 1
  for(I in LL[1]:UR[1]){
    for(J in LL[2]:UR[2]){
      tilePath = paste(path,paste(zoom,I,J,sep="/"),'.png',sep='')
      tileBounds = .tile2boundingBox(I,J,zoom)
	  tileExtent= .tile2extentMercator(I,J,zoom)
      tileData[[i]]=list(path=tilePath,bounds=tileBounds,I=I,J=J,zoom=zoom,src=path,
			  extent=tileExtent)
      i = i + 1
    }
  }
  return(tileData)
}

getTiles <- function(xlim,ylim,zoom,path,maxTiles = 16,cacheDir=tempdir(),
		timeOut=5*24*60*60,verbose=FALSE){
  nt = nTiles(xlim,ylim,zoom)
  if(nt > maxTiles){
    stop("Cant get ",nt," tiles with maxTiles set to ",maxTiles)
  }
  tileData = getTilePaths(xlim,ylim,zoom,path)
  rasters = list()
  localStore = FALSE
  if(file.exists(path)){
    if(!file.info(path)$isdir){
      stop("Path ",path," is not a folder")
    }
    localStore = TRUE
  }

  thecrs =CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")		


  colourtable = NULL
  for(ip in 1:length(tileData)){
    p = tileData[[ip]]$path
    if(localStore){
      where = p
    }else{
      where = .getTileCached(tileData[[ip]],cacheDir,timeOut=timeOut,
			  verbose=verbose)
    }
	
	thisimage = raster(where)
	newcolours = thisimage@legend@colortable
	newcolours = newcolours[!newcolours %in% names(colourtable)]
	toadd = seq(length(colourtable), len=length(newcolours),by=1)
	names(toadd) = newcolours
	colourtable = c(colourtable, toadd)
	newvalues = colourtable[thisimage@legend@colortable]
	
	values(thisimage) = newvalues[values(thisimage)+1]
	extent(thisimage) = tileData[[ip]]$extent
	proj4string(thisimage) = thecrs
	
	rasters[[ip]] = thisimage

	
	}

	rasters = do.call(merge, rasters)
	rasters@legend@colortable = names(colourtable)
	
	return(rasters)	
	}
	

.getTileCached <- function(tileData,cacheDir,timeOut=30,verbose=FALSE){
  srcDir = paste("X",make.names(tileData$src),sep="")
  tileDir = paste(file.path(cacheDir,srcDir,tileData$zoom,tileData$I))
  tilePath = paste(file.path(tileDir,tileData$J),".png",sep="")
  retrieve = FALSE
  if(!file.exists(tilePath)){
    if(verbose)cat("Getting new tile\n")
    retrieve = TRUE
  }else{
    age = difftime(Sys.time(),file.info(tilePath)$mtime,units = "mins")
    if (age > timeOut){
      if(verbose)cat("Tile aged ",age," expired from cache\n")
      retrieve = TRUE
    }else{
      if(verbose){cat("Tile found in cache\n")}
    }
  }
  if(retrieve){
    if(verbose)cat("Downloading tile\n")
    dir.create(tileDir,recursive=TRUE,showWarnings=FALSE)
    p = tileData$path
    download.file(p,tilePath,mode="wb",quiet=!verbose)
  }
  
  return(tilePath)
}
