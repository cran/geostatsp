
pcPriorRange = function(q, p=0.5, cellSize=1) {
  if(!is.numeric(cellSize)) {
    cellSize = xres(cellSize)
  } 
  # suppose lower 5% quantile of rangeKm is 50km
# ... of range is 50000/cellSize
# upper 5% quantile (lower 95%) of scale is  cellSize / 50000 
  # 1-exp(-lambda x) = 1-p
  # 1-0.05 = 1-exp(-lambda * cellSize/50000)
# log(0.05) = -lambda * cellSize / 50000
# lambda = -log(0.05) * 50000/cellSize
  # lambda is a multiplicative parameter
  # I call lambda a scale, wikipedia and R call it rate
  quantileCells = q/cellSize
  lambda = -log(p) * quantileCells 
  top99 = stats::qexp(0.99, lambda)
  xSeqScale = seq(0, top99, len=1001)
  bottom99 = stats::qexp(0.05, lambda)
  xSeqRangeCells = seq(0, 1/bottom99, len=1001)
  xSeqRange = xSeqRangeCells * cellSize
  result = list(
    lambda = lambda,
    priorScale = cbind(
      x = xSeqScale, y = stats::dexp(xSeqScale, lambda)
      )
    )
  result$priorRange = cbind(
    x = xSeqRange,
    y = stats::dexp(1/xSeqRangeCells, lambda) * (xSeqRange)^(-2) * cellSize
    )  
  result$inla = paste0(
      "expression:
        lambda = ", lambda, ";
        scale = exp(-log_range);
        logdens = -lambda*scale + log(lambda);
        log_jacobian = - log_range;
        return(logdens + log_jacobian);", sep='')
  result$string = paste0(
    "prior=\"",
    result$inla,
    "\""  
    )
  result
}