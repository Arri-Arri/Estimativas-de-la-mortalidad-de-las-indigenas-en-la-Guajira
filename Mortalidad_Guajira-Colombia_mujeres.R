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

N = c(
  10464,
  29483,
  44899,
  40871,
  37622,
  34104,
  29784,
  25717,
  23516,
  19307,
  16513,
  13366,
  11689,
  8302,
  6874,
  4993,
  3408,
  4221)




D = c(
  204,
  58,
  23,
  26,
  46,
  43,
  44,
  50,
  57,
  65,
  64,
  94,
  101,
  127,
  135,
  139,
  175,
  571)


mx <- D/N
log(mx)


#Tabla de vida modelo realizada por el DAne colomnboa Ya mejorada*
standm <- c(
  -4.01961,
  -7.13962,
  -7.95881,
  -7.82566,
  -7.51659,
  -7.27313,
  -7.09296,
  -6.92157,
  -6.69926,
  -6.41946,
  -6.07237,
  -5.66508,
  -5.21056,
  -4.72218,
  -4.21757,
  -3.70442,
  -3.18858,
  -2.22632)

#Para calcular la tasa de mortalidad corregida cual es el denomidador? A mitad de periodo o puedo tomar el supuesto de población a mitad de periodo 

#Tabla completa de Guajira
std = c(
  0.04823,
  0.00479,
  0.00199,
  0.00114,
  0.00076,
  0.00059,
  0.00052,
  0.00050,
  0.00050,
  0.00054,
  0.00059,
  0.00065,
  0.00071,
  0.00079,
  0.00085,
  0.00091,
  0.00097,
  0.00102,
  0.00108,
  0.00113,
  0.00118,
  0.00124,
  0.00129,
  0.00134,
  0.00139,
  0.00144,
  0.00149,
  0.00154,
  0.00160,
  0.00166,
  0.00172,
  0.00179,
  0.00185,
  0.00193,
  0.00200,
  0.00208,
  0.00216,
  0.00225,
  0.00236,
  0.00246,
  0.00259,
  0.00273,
  0.00288,
  0.00305,
  0.00324,
  0.00345,
  0.00369,
  0.00394,
  0.00424,
  0.00456,
  0.00492,
  0.00533,
  0.00577,
  0.00627,
  0.00681,
  0.00739,
  0.00803,
  0.00872,
  0.00948,
  0.01031,
  0.01123,
  0.01224,
  0.01335,
  0.01457,
  0.01590,
  0.01739,
  0.01900,
  0.02082,
  0.02282,
  0.02503,
  0.02744,
  0.03010,
  0.03302,
  0.03620,
  0.03974,
  0.04367,
  0.04804,
  0.05288,
  0.05825,
  0.06415,
  0.07061,
  0.07764,
  0.08527,
  0.09350,
  0.10228,
  0.11158,
  0.12135,
  0.13163,
  0.14232,
  0.15336,
  0.16467,
  0.17618,
  0.18778,
  0.19938,
  0.21088,
  0.22218,
  0.23317,
  0.24372,
  0.25375,
  0.26313)

std <- log(std) ##Tabla de vida de la Guajira##



names(N) = names(D) = head(boundaries,-1)



## single-year log mortality rates from HMD
## these are the targets for TOPALS estimation

#Tabla de vida de Colombia##
ITA_HMD_logmx =  c( 
  -4.0196,
  -6.5642,
  -7.1562,
  -7.5056,
  -7.7287,
  -7.8753,
  -7.9576,
  -7.9576,
  -8.0472,
  -7.9576,
  -7.9294,
  -7.9294,
  -7.8240,
  -7.7517,
  -7.7063,
  -7.6211,
  -7.5811,
  -7.5239,
  -7.4525,
  -7.4186,
  -7.3385,
  -7.3385,
  -7.2644,
  -7.2502,
  -7.1954,
  -7.1562,
  -7.1309,
  -7.0941,
  -7.0470,
  -7.0356,
  -6.9911,
  -6.9696,
  -6.9178,
  -6.8880,
  -6.8495,
  -6.7944,
  -6.7508,
  -6.7171,
  -6.6377,
  -6.6003,
  -6.5501,
  -6.4890,
  -6.4192,
  -6.3654,
  -6.2926,
  -6.2247,
  -6.1516,
  -6.0792,
  -6.0035,
  -5.9257,
  -5.8430,
  -5.7604,
  -5.6723,
  -5.5887,
  -5.4943,
  -5.4081,
  -5.3144,
  -5.2177,
  -5.1260,
  -5.0268,
  -4.9295,
  -4.8308,
  -4.7330,
  -4.6284,
  -4.5319,
  -4.4287,
  -4.3298,
  -4.2267,
  -4.1203,
  -4.0196,
  -3.9170,
  -3.8149,
  -3.7111,
  -3.6042,
  -3.5029,
  -3.3983,
  -3.2947,
  -3.1915,
  -3.0870,
  -2.9824,
  -2.8829,
  -2.7800,
  -2.6812,
  -2.5840,
  -2.4900,
  -2.3982,
  -2.3085,
  -2.2207,
  -2.1350,
  -2.0521,
  -1.9708,
  -1.8895,
  -1.8138,
  -1.7337,
  -1.6741,
  -1.6006,
  -1.5482,
  -1.4792,
  -1.4268,
  -1.3618)


show = function(fit, hue='#EF8737') {
  
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
         subtitle = paste(sum(fit$D),'deaths to',round(sum(fit$N)),'mujeres')) +
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

ggsave( filename  =  'C:/Users/Alí/Documents/PDh_Cedeplar/Segundo_year/Segundo semestre/Metodos indirectos/Grafica_2.jpg', width = 6, height = 4, dpi = 600)


age = c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)

plot(age, N, type='h', main = "Exposure by Age")

plot(age, log(mx), pch=16, main='True Log Mortality Rates by Age\nEstonia 2010')



B = bs(age, knots=c(0,1,10,20,40,70), degree=1 )


fitted_logmx = standm + B %*% fit$alpha

se_logmx = sqrt( diag (B %*% fit$covar %*% t(B)) )


plot( age, log(D/N), pch='+',cex=1.2,ylim=c(-10,0) ,
      main='La Guajira, Colombia, Mujeres [2018]\nTOPALS fit +95% CI')
rug(age[D==0], side=1, ticksize=.015)
lines(age, standm, lty=1, lwd=3, col='grey')
points(age, fitted_logmx, cex=.80, pch=16, col='#EF8737')
L = fitted_logmx - 1.96 * se_logmx
H = fitted_logmx + 1.96 * se_logmx
segments( age, L, age, H, col='#EF8737', lwd=.60)

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
         main='Simulación de los parametros \n Guajira mujeres, 2018')

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








