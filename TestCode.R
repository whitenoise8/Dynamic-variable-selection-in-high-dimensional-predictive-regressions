#: Last update 19 april 2023
#: Author: Nicolas Bianco (nicolas.bianco@phd.unipd.it)
#: [for any issue, feel free to mail me]

#: Example code for:
#: Dynamic variable selection in high-dimensional predictive regressions 
#: (Bernardi and Bianchi and Bianco, 2023)

#: Required packages: 
#: Rcpp, RcppArmadillo, RcppEigen, RcppNumerical, BH,
#: ggplot2, TeachingDemos, splines2
#: Download dependencies
list_of_packages = c("Rcpp", "RcppArmadillo", "RcppEigen", "RcppNumerical", "BH", 
                     "ggplot2", "TeachingDemos", "splines2")
new_packages = list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)>0) install.packages(new_packages)

#: Load package
library(dynBG)

# Libraries
library(splines2)

# Example
set.seed(190423)

n = 200
p1 = 1 #: Always included
p_ = 2 #: Dynamic sparsity
p0 = 1 #: Always zero
p = p1+p_+p0
eta2 = 0.1 

alpha = 0 #: Intercept
s2 = 0.25 #: Error variance

# Simulate data (in this code constant intercept and homoskedastic)
sim_data = sptvpsim(n,p1,p_,p0,alpha,s2,eta2,lambda=100)

y = sim_data$y
X = sim_data$X
beta = sim_data$beta
gamma = sim_data$gamma

#: Sparse TVP with smoothing
hyper = list(As=0.01,Bs=0.01, #: Prior sigma^2
             Ae=0.01,Be=0.01, #: Prior eta^2
             Ax=2.00,Bx=5.00, #: Prior xi^2
             S2=1.00,         #: Prior constant intercept
             k0=100)
lw = 10
dg = 3
W = bSpline(1:n, 
            knots=seq(lw,n-lw,by=lw), 
            degree=dg, 
            intercept=TRUE)
hyper$W = W
mod = BGTVP(y,X,hyper,
            icept="CI",              #: intercept CI=constant, TVI=time-varying
            options=list(sv=1,       #: sv (stochastic volatility) 1=yes, 0=no 
                         smooth=1),  #: smooth 1=yes, 0=n0
            Trace=1)                 #: Print progress

#: Plot the trajectories of regression coefficients
par(mfrow=c(2,2),mar=c(2,2,2,2))
for (j in 1:p) {
  plot(beta[j,],type='l',col="grey",lwd=2,main=paste0("beta_",j))
  points(mod$mu_q_gamma[j,]*mod$mu_q_beta[j,],type='l',col=2,lwd=2)
}
par(mfrow=c(1,1),mar=c(2,2,2,2))

#: Plot the trajectories of posterior inclusion probabilities
par(mfrow=c(2,2),mar=c(2,2,2,2))
for (j in 1:p) {
  plot(gamma[j,],type='l',col="grey",lwd=2,ylim=c(0,1),main=paste0("gamma_",j))
  points(mod$mu_q_gamma[j,],type='l',col=2,lwd=2)
}
par(mfrow=c(1,1),mar=c(2,2,2,2))

#: Plot the heatmap (x=time,y=beta_j) of regression coefficients
matrixplot(mod$mu_q_gamma*mod$mu_q_beta)
#: Plot the heatmap (x=time,y=gamma_j) of posterior inclusion probabilities
matrixplot(mod$mu_q_gamma)



