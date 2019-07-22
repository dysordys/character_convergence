
require(deSolve) ## for integrating ODE systems
require(tidyverse) ## for efficient data manipulation and plotting


## right hand sides of the dynamical equations (equations 11 and 12 from
## McPeek 2019 Am Nat; DOI: 10.1086/701629)
## Input:
## - time: the moment at which the right hand side is evaluated
##         (needed for compatibility; not used)
## - state: the state variables, coerced into a vector, in the order:
##          (R[1], R[2], N[1], N[2], zR[1], zR[2], zN[1], zN[2]),
##          where R[i] (N[i]) is the density of resource (consumer) i, and
##          zR[i] (zN[i]) is the trait mean of resource (consumer) i
## - pars: list of model parameters; see the definition of "params1" below
##         to see what entries it should have
## Output:
## - the vector of time derivatives of the dynamical variables, in the same
##   order as in the state vector (dR[1]/dt, dR[2]/dt, dN[1]/dt, ..., dzN[2]/dt)
eqs <- function(time, state, pars) {
    ## give short names to state variables (R[i], N[i], zR[i], zN[i])
    R <- state[1:2]; N <- state[3:4]; zR <- state[5:6]; zN <- state[7:8]
    ## initialize derivatives
    dRdt <- c(0, 0); dNdt <- c(0, 0); dzRdt <- c(0, 0); dzNdt <- c(0, 0)
    ## time derivative of resource 1's density (eq. 11 in McPeek 2019 Am Nat)
    dRdt[1] <- (pars$c0[1]*(1-pars$gamma[1]*((zR[1]-pars$zRc[1])^2+pars$PR[1]))-
                pars$d[1]*R[1]-
                sum(pars$a0*pars$beta*N*
                    exp(-(zN-zR[1])^2/(2*(pars$PR[1]+pars$PN+pars$beta^2)))/
                    sqrt(pars$PR[1]+pars$PN+pars$beta^2)))*R[1]
    ## time derivative of resource 2's density (eq. 11 in McPeek 2019 Am Nat)
    dRdt[2] <- (pars$c0[2]*(1-pars$gamma[2]*((zR[2]-pars$zRc[2])^2+pars$PR[2]))-
                pars$d[2]*R[2]-
                sum(pars$a0*pars$beta*N*
                    exp(-(zN-zR[2])^2/(2*(pars$PR[2]+pars$PN+pars$beta^2)))/
                    sqrt(pars$PR[2]+pars$PN+pars$beta^2)))*R[2]
    ## time derivative of consumer 1's density (eq. 11 in McPeek 2019 Am Nat)
    dNdt[1] <- (pars$b[1]*pars$a0[1]*pars$beta[1]*
                sum(R*exp(-(zN[1]-zR)^2/(2*(pars$PR+pars$PN[1]+pars$beta[1]^2)))/
                    sqrt(pars$PR+pars$PN[1]+pars$beta[1]^2))-
                pars$f0[1]*(1+pars$theta[1]*((zN[1]-pars$zNf[1])^2+pars$PN[1]))-
                pars$g[1]*N[1])*N[1]
    ## time derivative of consumer 2's density (eq. 11 in McPeek 2019 Am Nat)
    dNdt[2] <- (pars$b[2]*pars$a0[2]*pars$beta[2]*
                sum(R*exp(-(zN[2]-zR)^2/(2*(pars$PR+pars$PN[2]+pars$beta[2]^2)))/
                    sqrt(pars$PR+pars$PN[2]+pars$beta[2]^2))-
                pars$f0[2]*(1+pars$theta[2]*((zN[2]-pars$zNf[2])^2+pars$PN[2]))-
                pars$g[2]*N[2])*N[2]
    ## time derivative of resource 1's trait mean (eq. 12 in McPeek 2019 Am Nat)
    dzRdt[1] <- (-2*pars$c0[1]*pars$gamma[1]*(zR[1]-pars$zRc[1])-
                 sum(pars$a0*pars$beta*N*(zN-zR[1])*
                     exp(-(zN-zR[1])^2/(2*(pars$PR[1]+pars$PN+pars$beta^2)))/
                     (pars$PR[1]+pars$PN+pars$beta^2)^(3/2)))*pars$GR[1]
    ## time derivative of resource 2's trait mean (eq. 12 in McPeek 2019 Am Nat)
    dzRdt[2] <- (-2*pars$c0[2]*pars$gamma[2]*(zR[2]-pars$zRc[2])-
                 sum(pars$a0*pars$beta*N*(zN-zR[2])*
                     exp(-(zN-zR[2])^2/(2*(pars$PR[2]+pars$PN+pars$beta^2)))/
                     (pars$PR[2]+pars$PN+pars$beta^2)^(3/2)))*pars$GR[2]
    ## time derivative of consumer 1's trait mean (eq. 12 in McPeek 2019 Am Nat)
    dzNdt[1] <- (-pars$b[1]*pars$a0[1]*pars$beta[1]*
                 sum(R*(zN[1]-zR)*
                     exp(-(zN[1]-zR)^2/(2*(pars$PR+pars$PN[1]+pars$beta[1]^2)))/
                     (pars$PR+pars$PN[1]+pars$beta[1]^2)^(3/2))-
                 2*pars$f0[1]*pars$theta[1]*(zN[1]-pars$zNf[1]))*pars$GN[1]
    ## time derivative of consumer 2's trait mean (eq. 12 in McPeek 2019 Am Nat)
    dzNdt[2] <- (-pars$b[2]*pars$a0[2]*pars$beta[2]*
                 sum(R*(zN[2]-zR)*
                     exp(-(zN[2]-zR)^2/(2*(pars$PR+pars$PN[2]+pars$beta[2]^2)))/
                     (pars$PR+pars$PN[2]+pars$beta[2]^2)^(3/2))-
                 2*pars$f0[2]*pars$theta[2]*(zN[2]-pars$zNf[2]))*pars$GN[2]
    ## coerce derivatives into a vector and then a list, then return them
    return(list(c(dRdt, dNdt, dzRdt, dzNdt)))
}

## put results from integrating the ODEs in tidy form
organize_results <- function(sol) {
    dat <- sol %>% as.data.frame %>% as_tibble ## convert to tibble
    names(dat)[1] <- "time" ## first column is time
    names(dat)[2:5] <- paste0("density_", 1:4) ## then the four densities
    names(dat)[6:9] <- paste0("trait value_", 1:4) ## and the four trait means
    dat <- dat %>%
        gather("var", "v", 2:ncol(dat)) %>% ## normalize using key-value pairs
        ## separate the single column of data type (density or trait mean) and
        ## species number information into two columns
        separate(var, c("type", "species"), sep="_") %>%
        ## separate columns for densities and trait means for each species
        spread(type, v) %>%
        ## rename species from (1, 2, 3, 4) to (R1, R2, N1, N2)
        mutate(species=c("R1", "R2", "N1", "N2")[as.integer(species)]) %>%
        ## indicate trophic level in a separate column
        mutate(type=ifelse(species%in%c("R1", "R2"), "resource", "consumer")) %>%
        ## merge density and trait mean columns into key-value pairs
        gather(attribute, v, c("density", "trait value"))
    return(dat)
}

## parameters, coerced into list
params1 <- list(
    c0=c(3.0, 3.0),      ## maximum birth rate of resource species
    d=c(0.02, 0.02),     ## strength of resource intraspecific density dependence
    gamma=c(0.01, 0.01), ## rate of decrease of resource intrinsic growth function
    zRc=c(20.0, 30.0),   ## where resource intrinsic growth function is maximal
    zNf=c(25.0, 25.0),   ## where consumer intrinsic growth function is minimal
    GR=c(0.2, 0.2),      ## additive genetic variance of resource species
    ER=c(0.4, 0.4),      ## environmental variance of resource species
    GN=c(0.2, 0.2),      ## additive genetic variance of consumer species
    EN=c(30.0, 30.0),    ## environmental variance of resource species
    a0=c(0.5, 0.5),      ## maximum attack coefficient of consumer species
    b=c(0.1, 0.1),       ## consumer conversion efficiencies
    beta=c(5.0, 5.0),    ## width of Gaussian attack coefficient curve
    g=c(0.0, 0.0),       ## strength of consumer intraspecific density dependence
    theta=c(0.01, 0.01), ## rate of increase of consumer intrinsic growth function
    f0=c(1.2, 1.2)       ## minimum consumer intrinsic death rate
)
## total phenotypic variance = genetic + environmental
params1$PR <- params1$GR + params1$ER ## for resources
params1$PN <- params1$GN + params1$EN ## and for consumers
params1$Rinit <- c(40, 50) ## initial resource densities
params1$Ninit <- c(1, 10) ## initial consumer densities
params1$zRinit <- c(20, 30) ## initial resource trait means
params1$zNinit <- c(16, 34) ## initial consumer trait means
## initial conditions coerced into one vector
params1$ic <- c(params1$Rinit, params1$Ninit, params1$zRinit, params1$zNinit)
## two other parameter sets
params2 <- params1 ## parameter set 2: like set 1, except
params2$EN[2] <- 31.0 ## EN[2] is slightly perturbed
params2$PR <- params2$GR + params2$ER ## recalculate phenotypic variance (resource)
params2$PN <- params2$GN + params2$EN ## recalculate phenotypic variance (consumer)
params3 <- params1 ## parameter set 3: like set 1, except
params3$ic[7] <- 20 ## zNinit[1] is slightly perturbed


## time sampling points when solving ODEs
tmax <- 1600 ## time to integrate equations for
stepout <- tmax/500 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time

## solve ODEs and put results in tidy table
## First, solve with parameter set 1
sol1 <- ode(func=eqs, y=params1$ic, parms=params1, times=time) %>%
    organize_results %>%
    mutate(scenario="(A) default") %>% ## add column labeling parameter set
    mutate(s=ifelse(species=="R1", sqrt(params1$PR[1]), ## take square root of
             ifelse(species=="R2", sqrt(params1$PR[2]), ## trait variances, to get
             ifelse(species=="N1", sqrt(params1$PN[1]), ## standard deviations
                    sqrt(params1$PN[2])))))
## Second, solve with parameter set 2
sol2 <- ode(func=eqs, y=params2$ic, parms=params2, times=time) %>%
    organize_results %>%
    mutate(scenario="(B) perturbed trait variance") %>%
    mutate(s=ifelse(species=="R1", sqrt(params2$PR[1]),
             ifelse(species=="R2", sqrt(params2$PR[2]),
             ifelse(species=="N1", sqrt(params2$PN[1]),
                    sqrt(params2$PN[2])))))
## Third, solve with parameter set 3
sol3 <- ode(func=eqs, y=params3$ic, parms=params3, times=time) %>%
    organize_results %>%
    mutate(scenario="(C) perturbed initial conditions") %>%
    mutate(s=ifelse(species=="R1", sqrt(params3$PR[1]),
             ifelse(species=="R2", sqrt(params3$PR[2]),
             ifelse(species=="N1", sqrt(params3$PN[1]),
                    sqrt(params3$PN[2])))))
## merge solutions into one large tibble
sol <- bind_rows(sol1, sol2, sol3) %>%
    ## only show trait std dev s if the attribute in the row is trait value
    mutate(s=ifelse(attribute=="trait value", s, 0))


## plot results
sol %>%
    ggplot() +
    aes(x=time, y=v, ymin=v-s, ymax=v+s, colour=species, linetype=species,
        fill=species) +
    geom_line() +
    geom_ribbon(colour=NA, alpha=0.1) +
    facet_grid(attribute~scenario, scales="free_y", switch="y") +
    scale_colour_manual(name="species", values=c(1, 4, 6, 2)) +
    scale_fill_manual(name="species", values=c(1, 4, 6, 2)) +
    scale_linetype_manual(name="species", values=c(1, 1, 2, 2)) +
    theme_bw() + theme(axis.title.y=element_blank(), legend.position="none")

