if(FALSE) {
  # creating the precision matrix entries in maternGmrfPrec.R
  
  params=c(oneminusar=0.1, shape=0, conditionalVariance=1)
  
  Ngrid = 5
  myGrid = geostatsp::squareRaster(extent(-Ngrid,Ngrid,-Ngrid,Ngrid), 1+2*Ngrid)
  midcell = cellFromRowCol(myGrid, Ngrid+1, Ngrid+1)
  
  precMat = maternGmrfPrec(myGrid, param=params)
  
  
  params1 = params
  params1['shape'] = 1
  precMat1a = maternGmrfPrec(myGrid, param=params1)
  precMat1m = precMat %*% precMat

  table(as.matrix(precMat1a)[,midcell])
  table(as.matrix(precMat1m)[,midcell])
  
  precMatChar = as.character(precMat)
  precMatChar[precMatChar==as.character(-1/attributes(precMat)$info$theo['a'])] = '-a'
  precMatChar = matrix(precMatChar, nrow(precMat), ncol(precMat))
  
  precForYacas = apply(precMatChar,1,paste,collapse=',')
  precForYacas = paste(precForYacas, collapse='},{')
  yacas0 = Ryacas::yacas(paste("prec0:= { {", precForYacas, "} }"))

  # matern 0
  shape = 0
  yacas0s = Ryacas::yacas(paste("prec", shape, sep=''))
  yacasExp = yacas0s$text
  
  yacasChar = gsub("\n", "", as.character(yacasExp))
  
  exprMat = as.matrix(read.table(
          text=gsub("([)][,] |list[(])list[(]", '\n', as.character(yacasChar)), 
          stringsAsFactors=FALSE, sep=','))
  nnIndexMat = as.matrix(NNmat(myGrid, nearest = shape+1))
  precEntriesMat = data.frame(ind=as.vector(nnIndexMat[,midcell]), eqn=as.vector(exprMat[,midcell]))
  precEntriesMat = precEntriesMat[precEntriesMat$ind >0, ]
  precEntriesMat = precEntriesMat[!duplicated(precEntriesMat$ind), ]
  precEntriesMat = precEntriesMat[order(as.integer(precEntriesMat$ind)),]
  precEntriesMat$ind = paste("'", precEntriesMat$ind, "'=", sep='')
  precEntriesMat$last = ','
  precEntriesMat[nrow(precEntriesMat), 'last'] = ''
  
  cat(paste(apply(precEntriesMat, 1, paste, collapse=' '), collapse='\n'), '\n')
  
  
  # matern 1, or 2nd order GRF  
  shape = 1
  yacas1 = Ryacas::yacas(paste("prec", shape, ":=MatrixPower(prec0,", shape+1, ")", sep=''))
  yacas1s = Ryacas::yacas(paste("prec", shape, sep=''))
  yacasExp = yacas1s$text
  
  yacasChar = gsub("\n", "", as.character(yacasExp))
  
  exprMat = as.matrix(read.table(
          text=gsub("([)][,] |list[(])list[(]", '\n', as.character(yacasChar)), 
          stringsAsFactors=FALSE, sep=','))
  nnIndexMat = as.matrix(NNmat(myGrid, nearest = shape+1))
  precEntriesMat = data.frame(ind=as.vector(nnIndexMat[,midcell]), eqn=as.vector(exprMat[,midcell]))
  precEntriesMat = precEntriesMat[precEntriesMat$ind >0, ]
  precEntriesMat = precEntriesMat[!duplicated(precEntriesMat$ind), ]
  precEntriesMat = precEntriesMat[order(as.integer(precEntriesMat$ind)),]
  precEntriesMat$ind = paste("'", precEntriesMat$ind, "'=", sep='')
  precEntriesMat$last = ','
  precEntriesMat[nrow(precEntriesMat), 'last'] = ''
  
  cat(paste(apply(precEntriesMat, 1, paste, collapse=' '), collapse='\n'), '\n')
  
  
  # matern 2, or 3rd order GRF  
  shape = 2
  yacas2 = Ryacas::yacas(paste("prec", shape, ":=MatrixPower(prec0,", shape+1, ")", sep=''))
  yacas2s = Ryacas::yacas(paste("Simplify(prec", shape, ")", sep=''))
  yacasExp = yacas2s$text
  
  yacasChar = gsub("\n", "", as.character(yacasExp))
  
  exprMat = as.matrix(read.table(
          text=gsub("([)][,] |list[(])list[(]", '\n', as.character(yacasChar)), 
          stringsAsFactors=FALSE, sep=','))
  nnIndexMat = as.matrix(NNmat(myGrid, nearest = shape+1))
  precEntriesMat = data.frame(ind=as.vector(nnIndexMat[,midcell]), eqn=as.vector(exprMat[,midcell]))
  precEntriesMat = precEntriesMat[precEntriesMat$ind >0, ]
  precEntriesMat = precEntriesMat[!duplicated(precEntriesMat$ind), ]
  precEntriesMat = precEntriesMat[order(as.integer(precEntriesMat$ind)),]
  precEntriesMat$ind = paste("'", precEntriesMat$ind, "'=", sep='')
  precEntriesMat$last = ','
  precEntriesMat[nrow(precEntriesMat), 'last'] = ''
  
  cat(paste(apply(precEntriesMat, 1, paste, collapse=' '), collapse='\n'), '\n')
  
  # matern 4, or 5th order GRF  
  shape = 4
  yacas4 = Ryacas::yacas(paste("prec", shape, ":=MatrixPower(prec0,", shape+1, ")", sep=''))
  yacas4s = Ryacas::yacas(paste("Simplify(prec", shape, ")", sep=''))
  yacasExp = yacas4s$text
  
  yacasChar = gsub("\n", "", as.character(yacasExp))
  exprMat = as.matrix(read.table(
          text=gsub("([)][,] |list[(])list[(]", '\n', as.character(yacasChar)), 
          stringsAsFactors=FALSE, sep=','))
  nnIndexMat = as.matrix(NNmat(myGrid, nearest = shape+1))
  precEntriesMat = data.frame(ind=as.vector(nnIndexMat[,midcell]), eqn=as.vector(exprMat[,midcell]))
  precEntriesMat = precEntriesMat[precEntriesMat$ind >0, ]
  precEntriesMat$eqn = gsub("[[:space:]]+", " ", precEntriesMat$eqn)
  precEntriesMat = precEntriesMat[!duplicated(precEntriesMat$eqn), ]
  precEntriesMat = precEntriesMat[order(as.integer(precEntriesMat$ind)),]
  precEntriesMat$ind = paste("'", precEntriesMat$ind, "'=", sep='')
  precEntriesMat$last = ','
  precEntriesMat[nrow(precEntriesMat), 'last'] = ''
  
  cat(paste(apply(precEntriesMat, 1, paste, collapse=' '), collapse='\n'), '\n')

  
  # testing

  Ngrid = 30
  myGrid = geostatsp::squareRaster(extent(-Ngrid,Ngrid,-Ngrid,Ngrid), 1+2*Ngrid)
  midcell = cellFromRowCol(myGrid, Ngrid+1, Ngrid+1)
  
  params4 = params0 = c(conditionalVariance=1, oneminusar=0.1, shape=0)
#  params4 = params0 = c(variance=100, range = 3, shape=0)
  params4['shape'] = 4
  
  precMat = maternGmrfPrec(myGrid, param=params0)
  precMatA = maternGmrfPrec(N= myGrid, param=params4)
  
  precMat4m = precMat %*% precMat %*% precMat %*% precMat %*% precMat
  
  
  subA = matrix(as.matrix(precMatA)[,midcell], nrow(myGrid), ncol(myGrid))
  subM = matrix(as.matrix(precMat4m)[,midcell], nrow(myGrid), ncol(myGrid))  
  
  cSeq = 25:35
  subA[cSeq, cSeq]
  subM[cSeq, cSeq]
  quantile(  subA[cSeq, cSeq] -  subM[cSeq, cSeq])

 
  
  }