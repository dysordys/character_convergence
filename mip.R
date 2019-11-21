
require(deSolve)
require(tidyverse)


## right hand sides of the dynamical equations (equations 11 and 12 from
## McPeek 2019 Am Nat; DOI: 10.1086/701629)
## Input:
## - time: the moment at which the right hand side is evaluated
##         (needed for compatibility; not used)
## - state: the state variables, coerced into a vector, in the order:
##          (R[1], R[2], N[1], N[2], zR[1], zR[2], zN[1], zN[2]),
##          where R[i] (N[i]) is the density of resource (consumer) i, and
##          zR[i] (zN[i]) is the trait mean of resource (consumer) i
## - pars: list of model parameters
## Output:
## - the vector of time derivatives of the dynamical variables, in the same
##   order as in the state vector (dR[1]/dt, dR[2]/dt, dN[1]/dt, ..., dzN[2]/dt)
eqs <- function(time, state, pars) {
    R <- state[1:2]
    N <- state[3:4]
    dRdt <- c(0, 0)
    dNdt <- c(0, 0)
    dRdt[1] <-(pars$c0[1]*(1-pars$gamma[1]*(
        (pars$zR[1]-pars$zRc[1])^2+pars$PR[1]))-pars$d[1]*R[1]-
        sum(pars$a0*pars$beta*N*
            exp(-(pars$zN-pars$zR[1])^2/(2*(pars$PR[1]+pars$PN+pars$beta^2)))/
            sqrt(pars$PR[1]+pars$PN+pars$beta^2)))*R[1]
    dRdt[2] <- (pars$c0[2]*(1-pars$gamma[2]*(
        (pars$zR[2]-pars$zRc[2])^2+pars$PR[2]))-pars$d[2]*R[2]-
        sum(pars$a0*pars$beta*N*
            exp(-(pars$zN-pars$zR[2])^2/(2*(pars$PR[2]+pars$PN+pars$beta^2)))/
            sqrt(pars$PR[2]+pars$PN+pars$beta^2)))*R[2]
    dNdt[1] <- (pars$b[1]*pars$a0[1]*pars$beta[1]*
                sum(R*exp(-(pars$zN[1]-pars$zR)^2/
                          (2*(pars$PR+pars$PN[1]+pars$beta[1]^2)))/
                    sqrt(pars$PR+pars$PN[1]+pars$beta[1]^2))-
                pars$f0[1]*(1+pars$theta[1]*
                            ((pars$zN[1]-pars$zNf[1])^2+pars$PN[1]))-
                pars$g[1]*N[1])*N[1]
    dNdt[2] <- (pars$b[2]*pars$a0[2]*pars$beta[2]*
                sum(R*exp(-(pars$zN[2]-pars$zR)^2/
                          (2*(pars$PR+pars$PN[2]+pars$beta[2]^2)))/
                    sqrt(pars$PR+pars$PN[2]+pars$beta[2]^2))-
                pars$f0[2]*(1+pars$theta[2]*
                            ((pars$zN[2]-pars$zNf[2])^2+pars$PN[2]))-
                pars$g[2]*N[2])*N[2]
    return(list(c(dRdt, dNdt)))
}

## per capita growth rates of the two consumers
## Input:
## - state: the state variables, coerced into a vector, in the order:
##          (R[1], R[2], N[1], N[2], zR[1], zR[2], zN[1], zN[2]),
##          where R[i] (N[i]) is the density of resource (consumer) i, and
##          zR[i] (zN[i]) is the trait mean of resource (consumer) i
## - pars: list of model parameters
pgr <- function(state, pars) {
    R <- state[1:2]
    N <- state[3:4]
    r1 <- (pars$b[1]*pars$a0[1]*pars$beta[1]*
           sum(R*exp(-(pars$zN[1]-pars$zR)^2/
                     (2*(pars$PR+pars$PN[1]+pars$beta[1]^2)))/
               sqrt(pars$PR+pars$PN[1]+pars$beta[1]^2))-
           pars$f0[1]*(1+pars$theta[1]*
                       ((pars$zN[1]-pars$zNf[1])^2+pars$PN[1]))-
           pars$g[1]*N[1])
    r2 <- (pars$b[2]*pars$a0[2]*pars$beta[2]*
           sum(R*exp(-(pars$zN[2]-pars$zR)^2/
                     (2*(pars$PR+pars$PN[2]+pars$beta[2]^2)))/
               sqrt(pars$PR+pars$PN[2]+pars$beta[2]^2))-
           pars$f0[2]*(1+pars$theta[2]*
                       ((pars$zN[2]-pars$zNf[2])^2+pars$PN[2]))-
           pars$g[2]*N[2])
    return(c(r1, r2))
}

## parameters, coerced into list
params <- list(
    c0=c(3.0, 3.0),
    d=c(0.02, 0.02),
    gamma=c(0.01, 0.01),
    zRc=c(20.0, 30.0),
    zNf=c(25.0, 25.0),
    GR=c(0.2, 0.2),
    ER=c(0.4, 0.4),
    GN=c(0.2, 0.2),
    EN=c(30.0, 30.0),
    a0=c(0.5, 0.5),
    b=c(0.1, 0.1),
    beta=c(5.0, 5.0),
    f=c(0.15, 0.15),
    g=c(0.0, 0.0),
    theta=c(0.01, 0.01),
    f0=c(1.2, 1.2),
    zR=c(15.9, 34.1),
    zN=c(25.0, 25.0)
)
params$PR <- params$GR + params$ER
params$PN <- params$GN + params$EN

## time sampling points when solving ODEs
tmax <- 100 ## time to integrate equations for
stepout <- 0.1 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time

## mutual invasibility table
zaxis <- seq(19, 31, l=301)
invasion_table <- expand.grid(zN1=zaxis, zN2=zaxis) %>%
    as_tibble %>%
    mutate(rinv1=0, rinv2=0)

## obtain mutual invasibility values
percent <- 0
for (i in 1:nrow(invasion_table)) {
    params$zN <- as.numeric(invasion_table[i,1:2])
    ## can consumer 1 invade consumer 2?
    sol <- ode(func=eqs, y=c(48.9, 48.9, 0, 5), parms=params, times=time)
    invasion_table$rinv1[i] <- pgr(sol[nrow(sol), -1], params)[1]
    ## can consumer 2 invade consumer 1?
    sol <- ode(func=eqs, y=c(48.9, 48.9, 5, 0), parms=params, times=time)
    invasion_table$rinv2[i] <- pgr(sol[nrow(sol), -1], params)[2]
    ## progress indicator
    if (round(100*i/nrow(invasion_table))>percent) {
        percent <- round(100*i/nrow(invasion_table))
        cat(paste0(percent, "%", " \n"))
    }
}

## convert numerical invasion growth rates to verbal descriptors
invasion_table <- invasion_table %>%
    mutate(rinv1=ifelse(abs(rinv1)<1e-10, 0, rinv1),
           rinv2=ifelse(abs(rinv2)<1e-10, 0, rinv2)) %>%
    mutate(invasibility=ifelse((rinv1>0)&(rinv2>0), "1 & 2",
                        ifelse((rinv1>0)&(rinv2<0), "only 1",
                        ifelse((rinv1<0)&(rinv2>0), "only 2",
                        ifelse((rinv1<0)&(rinv2<0), "neither", "neutral"))))) %>%
    mutate(invasibility=factor(invasibility, ordered=TRUE,
                               levels=c("1 & 2", "only 1", "only 2", "neutral")))

## plot results
ggplot(invasion_table) +
    aes(x=zN1, y=zN2, fill=invasibility) +
    geom_raster(interpolate=TRUE) +
    scale_fill_manual(values=c("#377eb8", "gray77", "gray87", "gray94")) +
    scale_x_continuous("trait of species 1", expand=c(0, 0)) +
    scale_y_continuous("trait of species 2", expand=c(0, 0)) +
    theme_bw() +
    theme(panel.grid=element_blank())
