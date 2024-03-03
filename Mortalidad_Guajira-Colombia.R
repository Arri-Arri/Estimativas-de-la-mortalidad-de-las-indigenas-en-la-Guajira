#Estimacion de la mortalidad en indigenas de la Guajira##
##14-12-2023
##Ali Miguel Arrieta-Arrieta

rm(list=ls())
gc()
setwd("C:/Users/Alí/Documents/PDh_Cedeplar/Segundo_year/Segundo semestre/Metodos indirectos/Mort_Guajira")

install.packages("splines")

library(tidyverse, quietly = TRUE)
library(splines, quietly = TRUE)


TOPALS_fit = function( N, D, std,
                       age_group_bounds   = 0:100,
                       knot_positions     = c(0,1,10,20,40,70), 
                       penalty_precision  = 2,
                       max_iter           = 50,
                       alpha_tol          = .00005,
                       details            = FALSE) {
  
  require(splines)
  
  ## single years of age from 0 to (A-1)
  A   = length(std)
  age = 0:(A-1)
  
  ## B is an AxK matrix. Each column is a linear B-spline basis function
  B      = bs( age, knots=knot_positions, degree=1 )
  K = ncol(B) 
  
  D1 = diff( diag(K), diff=1)
  P  = penalty_precision * crossprod(D1)
  
  ## number and width of age groups
  G     = length(age_group_bounds)-1   
  nages = diff(age_group_bounds)
  
  ## weighting matrix for mortality rates (assumes uniform
  ## distribution of single-year ages within groups)
  W = matrix(0, nrow=G, ncol=A, 
             dimnames=list(head(age_group_bounds,-1) , age))
  
  offset = 0
  for (g in 1:G) {
    W[g, offset + 1:nages[g]] = 1/nages[g]
    offset = offset + nages[g]
  }
  
  ## penalized log lik function
  Q = function(alpha) {
    M = W %*% exp( std + B %*% alpha)
    likelihood = sum(D * log(M) - N * M)
    penalty    = 1/2 * t(alpha) %*% P %*% alpha
    return( likelihood - penalty )
  }
  
  #------------------------------------------------
  # iteration function: 
  # next alpha vector as a function of current alpha
  #------------------------------------------------
  next_alpha = function(alpha) {
    mu = as.vector( exp( std + B %*% alpha))
    M  = as.vector( W %*% mu)
    
    Dhat = N * M
    
    X = W %*% diag(mu) %*% B
    A = diag(N/M)
    
    y = (D-Dhat)/N + X %*% alpha
    
    updated_alpha = solve( t(X) %*% A %*% X + P, t(X) %*% A %*% y)
    return(as.vector(updated_alpha))
  }
  
  ## main iteration:     
  a = rep(0, K)
  
  niter = 0
  repeat {
    niter      = niter + 1
    last_param = a
    a          = next_alpha( a )  # update
    change     = a - last_param
    
    converge = all( abs(change) < alpha_tol)
    overrun  = (niter == max_iter)
    
    if (converge | overrun) { break }
    
  } # repeat
  
  if (details | !converge | overrun) {
    if (!converge) print('did not converge')
    if (overrun) print('exceeded maximum number of iterations')
    
    mu    = as.vector( exp(std + B %*% a))
    M     = as.vector( W %*% mu )
    dhat  = N * M
    
    X     = W %*% diag(mu) %*% B
    A     = diag(N/M)
    
    covar = solve( t(X) %*% A %*% X + P)
    
    return( list( alpha             = a, 
                  D                 = D,
                  N                 = N,
                  age_group_bounds  = age_group_bounds,
                  knots             = knot_positions,
                  std               = std,
                  B                 = B,
                  logm              = std + B %*% a,
                  covar             = covar,
                  Qvalue            = Q(a),
                  converge          = converge, 
                  maxiter           = overrun))
  } else return( a) 
  
} # TOPALS_fit
## Italian female 1980 HMD data for age groups


boundaries = c(0,1,seq(5,85,5))  # last group is [80,85)

N = c(10888,30335,47100,41875,37554,32886,27709,23706,21321,17687,15107,12222,10478,7584,6079,4239,2754,2955)

D = c(261.0,72.0,34,31,98,150,137,124,109,107,114,123,152,158,176,187,196,488)
mx <- D/N
log(mx)


#Tabla de vida modelo realizada por el DAne para Guajira* Ya mejorada*
standm <- c(
  -3.73012,
  -6.85919,
  -7.72728,
  -7.30393,
  -6.27573,
  -5.75002,
  -5.61723,
  -5.68750,
  -5.77782,
  -5.74246,
  -5.49615,
  -5.10046,
  -4.66352,
  -4.21560,
  -3.76471,
  -3.31790,
  -2.88145,
  -2.11233  )



#Para calcular la tasa de mortalidad corregida cual es el denomidador? A mitad de periodo o puedo tomar el supuesto de población a mitad de periodo 


std = c(0.06515,
0.00511,
0.00266,
0.00178,
0.00129,
0.00099,
0.00080,
0.00067,
0.00060,
0.00058,
0.00061,
0.00071,
0.00087,
0.00108,
0.00139,
0.00171,
0.00205,
0.00242,
0.00280,
0.00317,
0.00354,
0.00389,
0.00422,
0.00450,
0.00476,
0.00498,
0.00516,
0.00531,
0.00542,
0.00551,
0.00556,
0.00558,
0.00558,
0.00556,
0.00552,
0.00547,
0.00541,
0.00535,
0.00530,
0.00526,
0.00523,
0.00523,
0.00525,
0.00530,
0.00540,
0.00558,
0.00582,
0.00613,
0.00652,
0.00695,
0.00742,
0.00795,
0.00853,
0.00916,
0.00985,
0.01059,
0.01138,
0.01224,
0.01318,
0.01419,
0.01529,
0.01649,
0.01779,
0.01920,
0.02071,
0.02236,
0.02414,
0.02606,
0.02813,
0.03036,
0.03276,
0.03535,
0.03814,
0.04113,
0.04439,
0.04797,
0.05188,
0.05617,
0.06087,
0.06608,
0.07181,
0.07809,
0.08493,
0.09232,
0.10012,
0.10829,
0.11674,
0.12558,
0.13468,
0.14396,
0.15334,
0.16275,
0.17211,
0.18132,
0.19030,
0.19893,
0.20714,
0.21483,
0.22190,
0.22827)

std <- log(std) ##Tabla de vida de la Guajira##

std1 = c(-3.8933, -5.7776, -6.8474, -7.3298, -7.4519, -7.4408, -7.4807,
        -7.5845, -7.7219, -7.8628, -7.9771, -8.041, -8.0568, -8.0329,
        -7.9779, -7.9005, -7.8088, -7.7101, -7.6113, -7.5195, -7.4415,
        -7.3823, -7.3393, -7.308, -7.2837, -7.2619, -7.238, -7.2082,
        -7.1712, -7.1264, -7.0735, -7.0118, -6.9414, -6.8648, -6.7849,
        -6.7047, -6.6272, -6.5544, -6.4845, -6.4147, -6.3423, -6.2645,
        -6.1791, -6.0872, -5.9904, -5.8903, -5.7887, -5.6869, -5.586,
        -5.4866, -5.3895, -5.2953, -5.2049, -5.1186, -5.0347, -4.9513,
        -4.8664, -4.778, -4.6847, -4.5877, -4.4887, -4.3895, -4.2918,
        -4.1969, -4.1041, -4.0122, -3.9199, -3.8261, -3.7296, -3.6303,
        -3.5278, -3.4221, -3.3129, -3.2004, -3.0861, -2.9716, -2.8589,
        -2.7497, -2.6458, -2.5482, -2.4556, -2.3659, -2.2771, -2.187,
        -2.0942, -1.9991, -1.9028, -1.8062, -1.7105, -1.6164, -1.5242,
        -1.434, -1.3458, -1.2596, -1.1758, -1.0958, -1.0212, -0.9535,
        -0.8944, -0.8454)





aux <- 716/1423566.3
log(aux)

names(N) = names(D) = head(boundaries,-1)



## single-year log mortality rates from HMD
## these are the targets for TOPALS estimation

#Tabla de vida de Colombia##
ITA_HMD_logmx =  c( 
  -3.7301,
  -6.3368,
  -6.8685,
  -7.1691,
  -7.4021,
  -7.5617,
  -7.6843,
  -7.7753,
  -7.8240,
  -7.8240,
  -7.7287,
  -7.6009,
  -7.3858,
  -7.1691,
  -6.8782,
  -6.6687,
  -6.4567,
  -6.2818,
  -6.1193,
  -5.9955,
  -5.8889,
  -5.8091,
  -5.7353,
  -5.6840,
  -5.6493,
  -5.6213,
  -5.6103,
  -5.6103,
  -5.6130,
  -5.6324,
  -5.6380,
  -5.6665,
  -5.6869,
  -5.7138,
  -5.7322,
  -5.7572,
  -5.7667,
  -5.7828,
  -5.7893,
  -5.7926,
  -5.7828,
  -5.7699,
  -5.7509,
  -5.7260,
  -5.6810,
  -5.6352,
  -5.5780,
  -5.5041,
  -5.4307,
  -5.3538,
  -5.2726,
  -5.1904,
  -5.1077,
  -5.0237,
  -4.9351,
  -4.8460,
  -4.7607,
  -4.6692,
  -4.5795,
  -4.4936,
  -4.4022,
  -4.3095,
  -4.2220,
  -4.1320,
  -4.0393,
  -3.9492,
  -3.8594,
  -3.7679,
  -3.6773,
  -3.5885,
  -3.4963,
  -3.4082,
  -3.3204,
  -3.2289,
  -3.1426,
  -3.0517,
  -2.9670,
  -2.8790,
  -2.7925,
  -2.7059,
  -2.6188,
  -2.5334,
  -2.4468,
  -2.3644,
  -2.2818,
  -2.2030,
  -2.1272,
  -2.0546,
  -1.9841,
  -1.9168,
  -1.8489,
  -1.7872,
  -1.7263,
  -1.6676,
  -1.6142,
  -1.5672,
  -1.5127,
  -1.4761,
  -1.4223,
  -1.3900)


show = function(fit, hue='#41B7C4') {
  
  df_grouped = data.frame(
    L = head( fit$age_group_bounds, -1),
    U = tail( fit$age_group_bounds, -1),
    N = fit$N,
    D = fit$D
  ) %>%
    mutate(logmx_obs = log(D/N))
  
  
  df_single  = data.frame(
    age=  seq(fit$std) - .50,  # 0.5, 1.5, ...
    std = fit$std,
    logmx_true = ITA_HMD_logmx,
    logmx_fit  = fit$logm
  )
  
  this_plot =
    ggplot(data = df_single, aes(x=age,y=logmx_true)) +
    geom_line(aes(x=age,y=std), color='black', lwd=0.5) +
    geom_line(aes(x=age,y=logmx_fit), color=hue, lwd=3, alpha=.40) +
    geom_segment(data=df_grouped,aes(x=L,xend=U,
                                     y=logmx_obs,
                                     yend=logmx_obs),
                 color=hue,lwd=1.5, alpha=.90) +
    geom_point(size=0.60) +
    labs(x='Edad',y='Log (mx)',
         title='La Guajira, Colombia [2018]',
         subtitle = paste(sum(fit$D),'deaths to',round(sum(fit$N)),'hombres')) +
    scale_x_continuous(breaks=c(0,1,seq(5,100,5)),minor_breaks = NULL) +
    scale_y_continuous(limits=c(-10,0),breaks=seq(-10,0,2),minor_breaks = NULL) +
    theme_bw()
  
  print(this_plot)
} # show  



fit = TOPALS_fit(N,D,std,
                 age_group_bounds = boundaries,
                 details=TRUE)

str(fit)

show(fit)

ggsave( filename  =  'C:/Users/Alí/Documents/PDh_Cedeplar/Segundo_year/Segundo semestre/Metodos indirectos/Graphs/Grafica_1.jpg', width = 6, height = 4, dpi = 600)


age = c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)

plot(age, N, type='h', main = "Exposure by Age")

plot(age, log(mx), pch=16, main='True Log Mortality Rates by Age\nEstonia 2010')



B = bs(age, knots=c(0,1,10,20,40,70), degree=1 )


fitted_logmx = standm + B %*% fit$alpha

se_logmx = sqrt( diag (B %*% fit$covar %*% t(B)) )


plot( age, log(D/N), pch='+',cex=1.2,ylim=c(-10,0) ,
      main='La Guajira, Colombia [2018]\nTOPALS fit +95% CI')
rug(age[D==0], side=1, ticksize=.015)
lines(age, standm, lty=1, lwd=3, col='grey')
points(age, fitted_logmx, cex=.80, pch=16, col='#41B7C4')
L = fitted_logmx - 1.96 * se_logmx
H = fitted_logmx + 1.96 * se_logmx
segments( age, L, age, H, col='#41B7C4', lwd=.60)
ggsave( filename  =  'C:/Users/Alí/Documents/PDh_Cedeplar/Segundo_year/Segundo semestre/Metodos indirectos/Graphs/Grafica_2.jpg', width = 6, height = 4, dpi = 600)



##Simulacion de la mortalidad


CH = t(chol(fit$covar))
nsim = 1000
nshow = 20

sim_alpha = as.numeric(fit$alpha) +
            CH %*% matrix( rnorm(nsim*length(fit$alpha)),
                 ncol=nsim)

sim_lambda = standm + B %*% sim_alpha

matplot( age, sim_lambda[,1:nshow], type='l',ylim=c(-10,0),
         main='Simulación de los parametros \n Guajira hombres, 2018')

ggsave( filename  =  'C:/Users/Alí/Documents/PDh_Cedeplar/Segundo_year/Segundo semestre/Metodos indirectos/Grafica_3.jpg', width = 6, height = 4, dpi = 600)



#esperanza de vida 

e0 = function(logmx) {
  mx = exp(logmx)
  px = exp(-mx)
  lx = c(1,cumprod(px))
  return( sum(head(lx,-1) + tail(lx,-1)) / 2)
}

sim_e0 = apply(sim_lambda, 2, 'e0')

Q10 = quantile(sim_e0, .10)
Q50 = quantile(sim_e0, .50)
Q90 = quantile(sim_e0, .90)
round( head( sim_e0), 2)

## [1] 79.76 79.74 79.42 78.35 79.91 79.37
plot( density(sim_e0, adjust=1.5),
      main='Pará de Minas females 2010\nUncertainty about e0\nestimated from uncertainty about alpha')
points( Q50, .02, pch=16, cex=1.2)
segments(Q10, .02, Q90, .02, lwd=1.2)
text( Q50, .06, '80% interval')






e0 = function(logmx) {
  x = c(seq(logmx)-1, length(logmx)) # 0,1,...,start open interv
  mx = exp(c(logmx, tail(logmx,1))) # add open interval
  nx = c( rep(1, length(logmx)), Inf)
  px = exp(-mx * nx)
  lx = c(1,cumprod( head(px,-1))) # for 0,1,...,start open interv
  ax = 1/mx - 1*(px/(1-px))
  ax[1] = 0.1
  dx = -diff(c(lx,0))
  Lx = lx*px + dx*ax
  return( sum(Lx))
}


esim = apply(sim_lambda,2, e0)












##MUjeres

standm <- c(-1.31668, -2.66295, -3.27607, -3.14510, -2.99050, -2.89090, -2.81095, -2.73062, -2.64512, -2.53830, -2.40093, -2.23547,-2.05680,-1.87203,
            -1.67919,  -1.48022,  -1.27692, -0.91357)






