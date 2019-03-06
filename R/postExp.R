
postExp = function(x, 
  exclude = grep('^(range|aniso|shape|boxcox)', rownames(x)), 
  invLogit=FALSE) {

  if(is.numeric(exclude)) exclude = rownames(x)[exclude]
  
  if(is.logical(invLogit)) {
    if(invLogit) {
      invLogit = grep('[Ii]ntercept', rownames(x), value=TRUE)
    } else {
      invLogit = NULL
    }
  }
  if(is.numeric(invLogit)) invLogit = rownames(x)[invLogit]
  invLogit = setdiff(invLogit, exclude)

  colsIn = paste0(c(0.5,0.025, 0.975), 'quant')
  if(!all(colsIn %in% colnames(x))) {
    colsIn = c('estimate', 'ci0.025', 'ci0.975')
    if(!all(colsIn %in% colnames(x))) {
      warning('cant find quantiles in columns of x')
    }
  }

  colsOut = c('Estimate','2.5','97.5')

  if('Estimated' %in% colnames(x)) {
    x = x[x[,'Estimated'],]
  }

  res = x[,colsIn]
  colnames(res) = colsOut

  res[setdiff(rownames(res), exclude), ] = exp(
    res[setdiff(rownames(res), exclude), ])

  sdParams = grep('^sd([[:space:]]|$|Obs$|Nugget$|Spatial$)', rownames(res), value=TRUE)
  sdParams = setdiff(sdParams, exclude)

  rownames(res) = gsub(
    paste0('^(', paste(sdParams, collapse='|'), ')$'),
    'exp(\\1)', rownames(res)
  )
  if(length(sdParams) > 1)
    rownames(res) = gsub('exp[(]sd[)]', 'exp(sd spatial)', rownames(res))

  res[invLogit,] = res[invLogit,] / (1 + res[invLogit,])

  attributes(res)$caption =  'Posterior medians and credible intervals.' 
  attributes(res)$invLogit = invLogit
  attributes(res)$exclude = exclude

  res
}
