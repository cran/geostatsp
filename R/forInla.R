lcOneRow = function(thisrow, idxCol=NULL) {
  thisrow = thisrow[!is.na(thisrow)]
  if(length(thisrow)) {
    thisrow = sapply(thisrow, function(qq) list(list(weight=qq)))
    for(D  in idxCol)
      thisrow[[D]] = list(
        weight=1, 
        idx=thisrow[[D]]$weight
        )
    for(D in names(thisrow))
      thisrow[[D]] = thisrow[D]
    names(thisrow) = paste("v", 1:length(thisrow), sep="")
  }
  thisrow
}

inla.models=function(){
  if(requireNamespace("INLA", quietly=TRUE)){
    return(INLA::inla.models())
  } else {
    return(NULL)
  }
}