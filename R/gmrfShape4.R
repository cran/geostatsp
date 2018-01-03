if(FALSE) {
  # creating the precision matrix entries in maternGmrfPrec.R
  

  shape = Dnu = 5
  Dorder = Dnu + 1
  DnuEx = pmax(Dorder, 3)
  DD = geostatsp::NNmat(2*DnuEx+1, nearest=1)
  DD@x = c(-4,1)[DD@x]
  Ryacas::yacas(paste("DD := {{", paste(apply(DD,1,paste,collapse=','), collapse='},{'), "}}"))
  Ryacas::yacas("precCar1 := (1-theta) * Identity(Length(DD)) - (theta/4) * DD")
  Ryacas::yacas(paste(c(
						  "precCarNu := {BaseVector(Ceil(Length(precCar1)/2), Length(precCar1))}",rep('precCar1', Dorder)
				  ), collapse=' * '))
  Ryacas::yacas("precSubset := Partition(precCarNu[1], Sqrt(Length(precCar1)))")
  
  Ryacas::yacas("Echo(Simplify(precSubset Where theta == 4*a))")
  
  res = Ryacas::yacas("Echo(Simplify(precSubset Where theta == 4*a))")

  exprMat = as.matrix(read.table(
				  text=gsub('[{]|[}]', '', gsub("[}] [{]", '\n', as.character(res$PrettyForm))), 
				  stringsAsFactors=FALSE, sep=',', nrows = 2*DnuEx+1))
  
  nnIndexMat = geostatsp::NNmat(2*DnuEx+1, nearest = Dorder)
  nnIndexMat = matrix(nnIndexMat[,ceiling(ncol(nnIndexMat)/2)], 2*DnuEx+1) 
  nnIndexMat	  
  
  res =data.frame(as.vector(nnIndexMat), as.vector(exprMat))
  res = res[!duplicated(res[,1]) & res[,1] != '0',]
  res = res[order(res[,1]),]

cat(paste('\n"', res[,1], '" = ', res[,2], ',',sep=''))
  

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