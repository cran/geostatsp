nnIndices = function(index, baseDist) {
  rotDist =  Mod(baseDist)*exp(1i*(Arg(baseDist)+seq(0, 2*pi, len=5)))
  rotDist = c(rotDist, Im(rotDist)+1i*Re(rotDist))
  rotDist = unique(round(rotDist))
  
  cbind(index, x=Re(rotDist), y=Im(rotDist))
}


nbMat = rbind(
    c(1,0,0),# self
    c(2,0,1), #up
    c(2,0,-1), #down
    c(2,1,0), #right,
    c(2,-1,0), #left
    c(3,-1,-1),
    c(3,-1,1),
    c(3,1,-1),
    c(3,1,1),
    c(4,0,2), #up
    c(4,0,-2), #down
    c(4,2,0), #right,
    c(4,-2,0), #left
    c(5,-1,-2),
    c(5,-1,2),
    c(5,1,-2),
    c(5,1,2),
    c(5,-2,-1),
    c(5,-2,1),
    c(5,2,-1),
    c(5,2,1),
    c(6,0,3), #up
    c(6,0,-3), #down
    c(6,3,0), #right,
    c(6,-3,0), #left
    c(7,-2, 2), 
    c(7,-2,-2), 
    c(7, 2, 2), 
    c(7, 2,-2), 
    nnIndices(8,baseDist = 3+1i),
    nnIndices(9,baseDist = 4+0i),
    nnIndices(10,baseDist = 3+2i),
    nnIndices(11,baseDist = 4+1i),
    nnIndices(12,baseDist = 5+0i),
nnIndices(13,baseDist = 3+3i),
nnIndices(14,baseDist = 4+2i),
nnIndices(15,baseDist = 5+1i),
nnIndices(16,baseDist = 6+0i)
)


NNmat = function(N,Ny=N, nearest=3, adjustEdges=FALSE) {
  UseMethod("NNmat")	
}

NNmat.SpatRaster = function(N, Ny=N, nearest=3, adjustEdges=FALSE) {
  res = NNmat(ncol(N),nrow(N), nearest, adjustEdges)
  
  attributes(res)$raster= rast(N)
  
  res
}	

NNmat.default = function(N, Ny=N, nearest=3, adjustEdges=FALSE) {
  
  if(nearest<4){
    nbMat = nbMat[nbMat[,1] <= (c(2,4,6)[nearest]),] 
  }
  
  Nx = N
  theraster = rast(extent = ext(0,Nx, 0, Ny), nrows = Ny, ncol = Nx)
  Ncell = ncell(theraster)
  cellSeq = 1:Ncell
  terra::values(theraster) = cellSeq
  # find id's of cells on the border (for edge correction)  

 if(any(dim(theraster)[1:2] <= 2*nearest)) {
  stop("grid too small for this many neighbours") 
 }
  innerCells = crop(theraster, 
    ext(as.vector(ext(theraster)) + nearest * c(1,-1,1,-1))
  )
  
  innerCells = sort(values(innerCells))
  outerCells = sort(setdiff(values(theraster), innerCells))
  
  
  xMat = colFromCell(theraster, cellSeq)
  yMat = rowFromCell(theraster, cellSeq)
  
  xNb = outer(nbMat[,2],xMat, FUN='+')   
  yNb = outer(nbMat[,3],yMat, FUN='+')
  
  # which neighbour pairs are outside of the grid
  nbOutsideRegion = which(Matrix(xNb < 1) | Matrix(xNb > Nx) | 
          Matrix(yNb < 1) | Matrix(yNb> Ny))
  
  xNb = xNb[-nbOutsideRegion]
  yNb = yNb[-nbOutsideRegion]
  
  nbPairs = cbind(
      i = xNb + Nx*(yNb-1),
      j = matrix(cellSeq, nrow(nbMat), Ncell, byrow=TRUE)[-nbOutsideRegion],
      x = matrix(nbMat[,1], nrow(nbMat), Ncell)[-nbOutsideRegion])
  
  if(adjustEdges) {
    # set entries requiring edge correction to missing
    # allocating the memory here instead of in maternGmrfPrec is faster
    outerCellsPairs = expand.grid(i=outerCells, j=outerCells, x=NA)
    nbPairs = rbind(
        nbPairs, outerCellsPairs
    )
  }
  
  nbPairs = nbPairs[order(nbPairs[,'j'], nbPairs[,'i']),]
  
  result = Matrix::sparseMatrix(
      i=as.integer(nbPairs[,'i']),
      j=as.integer(nbPairs[,'j']),
      x=as.numeric(nbPairs[,'x']),
      dims = c(Ncell, Ncell),
      dimnames = list(
          cellSeq, cellSeq
      ), 
      symmetric=FALSE
  )
  
  
  attributes(result)$Nx = Nx
  attributes(result)$Ny = Ny
  attributes(result)$adjustEdges = adjustEdges
  attributes(result)$nearest = nearest
  attributes(result)$raster = theraster
  attributes(result)$cells = list(inner=innerCells, outer=outerCells)
  
  return(result)
}