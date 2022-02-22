# LIBRARIES 
## Manipulate data
library(dplyr)
library(lubridate)
## Solve equations
library(deSolve)
## Plotting
library(ggplot2)


# FUNCTION FOR THE FOOD WEB MODEL 
atn <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    
    N <- state[1]
    Xi <- state[-1]
    
    #### FORCING FUNCTIONS #################################################
    
    # Nutrients
    Nflux <- ninfun(times) #- n.out
    Nlim <- N/(N + kN)
    # Carrying capacity
    # approx_Ks <- approxfun(Ks, rule = 2)
    # Ks <- approx_Ks(times)
    # Ks.t <- Ks
    #Klim <- Gi(state, cij, Ks.t, prod) 
    Klim <- 1
    # Light
    #approx_light <- approxfun(I, rule = 2)
    I.t <- approx_light(times)
    # Temperature
    #approx_temp <- approxfun(T.i, rule = 2)
    T.i <- approx_temp(times)
    
    #temperature
    c.t = Q10^((T.i - T0)/10)
    x.t = x.i * c.t
    yij.t = yij #* c.t
    cp.t <- temp_lim(T.i, Topt = T0.phyt, Tmax = maxt, Tsteepness = tsteep)
    
    h <- approx_therm(times)
    # light lim
    prod2 <- prod
    prod2[1] <- 1
    sumphyt <- sum(Xi * prod2)
    I.z <- I.t * exp(-(a.phy * (sumphyt/1e9) + a.bg*h))
    c.l <- I.z/(I0 + I.z)
    r.l <- r * apply(cbind(Nlim, c.l), 1, min) * c(cp.t)
    #r.l <- r * apply(cbind(Nlim, c.l, cp.t), 1, min)
    
    #### DYNAMICS  #################################################
    
    w2 <- apply(w, 2, function(x){if(sum(x) > 0){x/sum(x)}else{rep(0, length(x))}})
    
    Nuptake <- sum(r.l * Xi * Klim * (1 - s) * n.frac)
    dN <- Nflux + Xi[1] * mineralization - Nuptake - N * n.out
    #dynamics
    dB <- r.l * Xi * Klim * (1 - s) +                                   
      # producer growth
      fa * x.t * Xi * colSums((yij.t * Fij(B0, q, d, p, Xi, w2))) -                           
      # gain from eating
      rowSums(t(apply(yij.t * Fij(B0, q, d, p, Xi, w2), 1, 
                      function(row){row * x.t * Xi}))/eij, na.rm = T) - 
      # loss from being eaten
      (fm * x.t * Xi) -                                                                        
      # maintenance loss
      m * Xi +                                                                                   
      # mortality (fish only)   
      x.in
    
    # increases in DOC
    DOCdelay <- 1#(1-exp(-.45 * ceiling(times)))^3500
    dB[1] <- (dB[1] + (sum(t(apply(yij.t * Fij(B0, q, d, p, Xi, w2), 1,
                                   function(row){row * x.t * Xi}))/eij  * (1-eij), na.rm = TRUE) +
                         # egestion
                         max(sum(r.l * Xi * Klim * s), 0))) * (1 - DOCloss) - (Xi[1] * mineralization)
    
    
    #chla <- (1/(6.4+45*Gi(state, cij, Ks.t, prod)*(1-c(0,c.l))))*prod
    
    return(list(c(dN, dB), L = I.z, tem = T.i))#, N.lim = Nlim, light = c.l[2], tlim = cp.t[2]))
  })
}

# FUNCTION TO DESCRIBE TEMPERATURE LIMITATION IN PHYTOPLANKTON
temp_lim <- function(temp, Topt = 25, Tmax = 40, Tsteepness = 150){
  temp_opt <- exp(-(((Topt - temp)^2)/Tsteepness))
  penalty <- 1 - (1/(1 + exp(-temp + Tmax)))
  
  return(temp_opt * penalty)
}

# FUNCTION FOR THE FUNCTIONAL RESPONSE
Fij <- function(B0, q, d, p, B, w){
  numer <-  w * B^q
  halfsat <- B0^q
  interfere <- 0
  #sapply(1:length(B), function(k){sum(d[k,] * p[,k] * B[k] * B0[k,]^q)}) # d * p * B * B0^q
  res <- colSums(w * B^q)
  
  fres <- numer/t(apply(halfsat, 1, function(x) x + interfere + res))
  fres[is.nan(fres)] <- 0
  return(fres)
}

# FUNCTION TO SET EXTINCTION THRESHOLD
goExtinct <- function(times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0
    return(c(states))
  })
}

# FUNCTION TO RUN THE SIMULATION 
run_NPZ <- function(state, tmax, params, temperature, light, nutrient, plot = TRUE){
  t.strt <- Sys.time()
  params$ninfun <- approxfun(nutrient, rule = 2)
  params$approx_light <- approxfun(light, rule = 2)
  params$approx_temp <- approxfun(temperature, rule = 2)
  
  out_sNPZ <- ode(state, times = 1:tmax, func = atn, parms = params, 
                  events = list(func = goExtinct, time = 1:tmax))
  
  p <- reshape2::melt(out_sNPZ[,2:ncol(out_sNPZ)]) %>% 
    ggplot(aes(x = Var1 + ymd("2015-05-02"), value)) +  
    geom_line(color = "blue", size = 1.2) +
    facet_wrap(~Var2, scales = "free_y") + 
    labs(x = "Date", y = "Concentration (ug/L)") + theme_bw()
  
  if(plot){print(p)}
  
  t.end <- Sys.time()
  t.diff <- t.end - t.strt
  cat("The simulation took", t.diff, units(t.diff))
  return(out_sNPZ)
}
