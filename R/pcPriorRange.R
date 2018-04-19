
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
  lambda = -log(p) * quantileCells / log(2)
  resultInla = paste0(
    "expression:
    lambda = ", lambda, ";
    scale = exp(-log_range);
    logdens = -lambda*scale + log(lambda);
    log_jacobian = - log_range;
    return(logdens + log_jacobian);", sep='')

  result= list(
    string = paste0("list(prior=\"", resultInla, "\")" ),
    param = c(q = q, prob = p), 
    extra = list(
      lambda=unname(lambda),
      medianRange = unname(cellSize*lambda), 
      medianScale = unname(1/(cellSize*lambda)),
      meanScale = unname(1/(log(2)*cellSize*lambda))),
    info = 'exponential prior for scale (pc prior)',
    cellSize = unname(cellSize))

  result$dprior = list(
    range = eval(parse(text=paste0(
      'function(x) x^(-2)*stats::dexp(1/x, rate=', 
      1/result$extra$meanScale, ')'))),
    scale = eval(parse(text=paste0(
      'function(x) stats::dexp(x, rate=', 1/result$extra$meanScale, ')')))
    )
  environment(result$dprior$range) = baseenv()
  environment(result$dprior$scale) = baseenv()
  result
}

gammaSd = function(param) {

  param = addNames(param)
# sd is gamma
# log(sd) is loggamma
# log(prec) = -2*log(sd) is dist'n below

# a=1/3;b=3/4;stuff = rgamma(10000, a,b)
# hist(log(stuff^(-2)), prob=TRUE)
# xseq = seq(-20, 20, len=1001) 
# lines(xseq, 
#  exp(a*log(b)  - lgamma(a) -log(2) - 
#  (a-1)*(xseq/2)-b*exp(-xseq/2)-xseq/2))
  #                    f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)              

  if(FALSE) {
    bob2 <- function(log_precision) {

      a = param['shape']
      b = param['rate']

      dens1 = stats::dgamma(exp(-log_precision/2),
        shape = param['shape'], rate = param['rate'], log=TRUE) - 
        log_precision/2 - log(2)

    logDens = -lgamma(a) - log(2) + a*log(b) - ((a-1)/2)  * log_precision - 
      exp(-log_precision/2)*b - log_precision/2
      exp(cbind(dens1, logDens))
    }
      # sample sd
    mySample = stats::rgamma(10001,
        shape = param['shape'], 
        rate = param['rate'])

      # histogram of log range
    hist(log(mySample^(-2)), breaks=100, prob=TRUE)
    logSeq = seq(-5,25,len=1001)
    matlines(logSeq, bob2(logSeq), col=c('blue','red'))
  }

  result = list(string = paste0(
    "list(prior=\"expression:",
     "\n    am1d2 = ", (param['shape']-1)/2, 
    ";\n    b = ", param['rate'], 
    ";\n    lngammalambda =", 
    -lgamma(param['shape']) - log(2) - param['shape'] * log(param['rate']),
    ";\n    return(lngammalambda - am1d2 * log_precision -", 
    "b*exp(-log_precision / 2.0) - log_precision/2.0);\")"
    ),
  param = param,
  info = 'gamma prior for sd',
  dprior = eval(parse(text=paste0('function(x) stats::dgamma(x, shape=',
    param['shape'],',rate=', param['rate'], ')')))
  )
  environment(result$dprior) = baseenv()
  result

}

addNames = function(param) {
  if(is.null(names(param))) 
    names(param) = c('shape','rate')
  param
}

gammaScale = function(param, cellSize) {

  param = addNames(param)

    # log gamma dens = - a * log(s) - loggamma(a) +
      # (a-1) * log(x) -x/s

    # log gamma dens = a * log(r) - loggamma(a) +
      # (a-1) * log(x) -x*r

    # y = log scale

    # log gamma dens = a * log(r) - loggamma(a) +
      # (a-1) * y -exp(y)*r

    # z = log (1/scale)

    # log gamma dens = a * log(r) - loggamma(a) +
      # (a-1) / z  - r * exp(-z)

# gsl_sf_lngamma(xx);
  rangePrior = list(
    param = param,
    info = 'gamma prior for scale',
    string = paste0(
      "list(prior=\"expression:\n    lambda = ", param['shape'], 
      ";\n    scale = ", param['rate'], 
      ";\n    lngammalambda = - lgamma(lambda)",
      ";\n    return(-lngammalambda + lambda * log(scale) - (lambda - 1.0) * log_range - scale * exp(-log_range) - log_range);",
      "\")"),
    extra = list(
      userParam = param[c('shape','rate')] * c(1,cellSize)
      )
    )

  rangePrior$dprior = list(
    range = eval(parse(text=paste0(
      'function(x) x^(-2)*stats::dgamma(1/x, shape=',
      rangePrior$extra$userParam['shape'],
      ',rate=', 
      rangePrior$extra$userParam['rate'], ')'))),
    scale = eval(parse(text=paste0(
      'function(x) stats::dgamma(x, shape=',
      rangePrior$extra$userParam['shape'],
      ',rate=', 
      rangePrior$extra$userParam['rate'], ')')))
    )


# should add - lgamma(a) 

  if(FALSE) {
    bob2 <- function(log_range) {

      a = param['shape']
      r = param['rate']

      dens1 = stats::dgamma(exp(-log_range), 
        shape = param['shape'], rate = param['rate'], log=TRUE) - 
      log_range

      logDens = a * log(r) - lgamma(a) + 
      (a - 1.0) / log_range - r*exp(-log_range) - log_range

      exp(cbind(dens1, logDens))
    }
      # sampe scale
    mySample = rgamma(10000, shape = param['shape'], rate = param['rate'])
      # histogram of log range
    hist(log(1/mySample), breaks=100, prob=TRUE)
    logRangeSeq = seq(-5,5,len=1001)
    matlines(logRangeSeq, bob2(logRangeSeq), col=c('blue','red'))
  }

  rangePrior
}