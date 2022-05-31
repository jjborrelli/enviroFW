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
    
    # current salt
    salt.i <- approx_salt(times)
    
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
      x.in -
      # immigration in
      (Xi * fsalt(salt.i, half = ksal, q = salsens))
      # salt mortality 
    
    # increases in DOC
    DOCdelay <- 1#(1-exp(-.45 * ceiling(times)))^3500
    dB[1] <- (dB[1] + (sum(t(apply(yij.t * Fij(B0, q, d, p, Xi, w2), 1,
                                   function(row){row * x.t * Xi}))/eij  * (1-eij), na.rm = TRUE) +
                         # egestion
                         max(sum(r.l * Xi * Klim * s), 0))) * (1 - DOCloss) - (Xi[1] * mineralization)
    
    
    #chla <- (1/(6.4+45*Gi(state, cij, Ks.t, prod)*(1-c(0,c.l))))*prod
    
    return(list(c(dN, dB)))#, N.lim = Nlim, light = c.l[2], tlim = cp.t[2]))
  })
}

atn_2N <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    
    N <- state[1:2]
    Xi <- state[-c(1:2)]
    
    #### FORCING FUNCTIONS #################################################
    
    # Nutrients
    Nflux <- c(ninfun1(times), ninfun2(times)) #- n.out
    
    # where N is a vector of nutrients
    # kN is a two row matrix with each row i having kN for nutrient i and columns are species
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
    #prod2[1] <- 1
    sumphyt <- sum(Xi * prod2)
    I.z <- I.t * exp(-(a.phy * (sumphyt/1e9) + a.bg*h))
    c.l <- I.z/(I0 + I.z)
    r.l <- r * apply(rbind(Nlim, c.l), 2, min) * c(cp.t)
    
    # current salt
    salt.i <- approx_salt(times)
    
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
      x.in -
      # immigration in
      (Xi * fsalt(salt.i, half = ksal, q = salsens))
      # salt mortality 
      
    
    # increases in DOC
    DOCdelay <- 1#(1-exp(-.45 * ceiling(times)))^3500
    dB[1] <- (dB[1] + (sum(t(apply(yij.t * Fij(B0, q, d, p, Xi, w2), 1,
                                   function(row){row * x.t * Xi}))/eij  * (1-eij), na.rm = TRUE) +
                         # egestion
                         max(sum(r.l * Xi * Klim * s), 0))) * (1 - DOCloss) - (Xi[1] * mineralization)
    
    
    #chla <- (1/(6.4+45*Gi(state, cij, Ks.t, prod)*(1-c(0,c.l))))*prod
    
    return(list(c(dN, dB)))
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

fsalt <- function(x, half, q){
  # using q=1 for high sensitivity, 1.5 for med, 2 for low
  res <- x^q / ((x^q + half^q))
  res[q == 0] <- 0
  return(res)
}

# FUNCTION TO SET EXTINCTION THRESHOLD
goExtinct <- function(times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0
    return(c(states))
  })
}

# FUNCTION TO RUN THE SIMULATION 
run_NPZ <- function(state, tmax, params, temperature, light, nutrient, plot = TRUE, plot_clip = 1){
  t.strt <- Sys.time()
  params$ninfun <- approxfun(nutrient, rule = 2)
  params$approx_light <- approxfun(light, rule = 2)
  params$approx_temp <- approxfun(temperature, rule = 2)
  
  out_sNPZ <- ode(state, times = 1:tmax, func = atn, parms = params, 
                  events = list(func = goExtinct, time = 1:tmax))
  
  p <- reshape2::melt(out_sNPZ[plot_clip:nrow(out_sNPZ),2:ncol(out_sNPZ)]) %>% 
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

get_ppar <- function(mat){
  p <- mat * 0
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      p[i,j] <- sum(mat[,i] == 1 & mat[,j] == 1)/sum(mat[,i] == 1 | mat[,j] == 1)
    }
  }
  p[is.nan(p)] <- 0
  return(p)
}

spnam <- c("DOC", "Nanno", "Blue", "Brown", "Green", "Red", "Flag", "Rot", "Clad", "Cal", "Cyc", "Fsh")


adj_npz <- matrix(c(0,0,1,
                    0,0,1,
                    0,0,0), nrow = 3, byrow = TRUE)

colnames(adj_npz) <- c("DOC", "Phyto", "Zoop")
rownames(adj_npz) <- c("DOC", "Phyto", "Zoop")



adj_gnpz <- matrix(c(0,0,0,0,0,0,1,1,1,1,0,
                     0,0,0,0,0,0,1,1,1,0,0,
                     0,0,0,0,0,0,0,1,1,1,1,
                     0,0,0,0,0,0,0,1,1,1,1,
                     0,0,0,0,0,0,0,1,1,1,1,
                     0,0,0,0,0,0,0,1,1,1,1,
                     0,0,0,0,0,0,0,0,1,1,1,
                     0,0,0,0,0,0,0,0,0,0,1,
                     0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0), nrow = 11, byrow = TRUE)

colnames(adj_gnpz) <- spnam[-length(spnam)]
rownames(adj_gnpz) <- spnam[-length(spnam)]

adj_gnpzf <- matrix(c(0,0,0,0,0,0,1,1,1,1,0,0,
                      0,0,0,0,0,0,1,1,1,0,0,0,
                      0,0,0,0,0,0,0,1,1,1,1,0,
                      0,0,0,0,0,0,0,1,1,1,1,0,
                      0,0,0,0,0,0,0,1,1,1,1,0,
                      0,0,0,0,0,0,0,1,1,1,1,0,
                      0,0,0,0,0,0,0,0,1,1,1,0,
                      0,0,0,0,0,0,0,0,0,0,1,0,
                      0,0,0,0,0,0,0,0,0,0,0,1,
                      0,0,0,0,0,0,0,0,0,0,0,1,
                      0,0,0,0,0,0,0,0,0,0,0,1,
                      0,0,0,0,0,0,0,0,0,0,0,0), nrow = 12, byrow = TRUE)

colnames(adj_gnpzf) <- spnam
rownames(adj_gnpzf) <- spnam

#mlist <- readRDS("data/smWebMats.rds")

adj <- c(list(adj_npz, adj_gnpz, adj_gnpzf))#, mlist)

#plist <- readRDS("data/smWebPARAMs.rds")


pars_npz <- list(
  n.in = .3,
  ninfun1 = approxfun(data.frame(x = 1:10000, y = 0.5)),
  ninfun2 = approxfun(data.frame(x = 1:10000, y = 1)),
  n.out = 0.2,   
  x.in = c(0, 0, 0),
  kN = rbind(rep(5,3),rep(10,3)),
  n.frac = rbind(rep(0.009,3), rep(0.15,3)),
  mineralization = 0.00,
  prod = c(0,1,0),
  r = c(0, 1, 0),
  x.i = c(0,0,0.08),
  s = 0.2,
  q = 1.2, 
  Ks = 1, 
  fa = 0.4,
  fm = 0.1,
  yij = adj_npz*10,
  B0 = adj_npz * 1500,
  eij = adj_npz * 0.66,
  d = adj_npz*rbeta(length(adj_npz), .5,4),
  w = apply(adj_npz, 2, function(x){if(sum(x) > 0){x/sum(x)}else{rep(0, length(x))}}),
  p = get_ppar(adj_npz),
  m = c(0,0,0.05), 
  cij = diag(rep(1.8,3)),
  h = 10,
  approx_therm = approxfun(data.frame(x = c(1, 10000), y = c(10, 10)), rule = 2),
  a.phy = 1e5,
  a.bg = 0.35, 
  I0 = c(1, 20, 1),
  approx_light = approxfun(data.frame(x = 1:10000, y = 150), rule = 2),
  Q10 = 3,
  T0 = 12,
  T0.phyt = c(20, 25, 20),
  tsteep = c(100, 150, 100),
  maxt = c(40, 40, 40),
  approx_temp = approxfun(data.frame(x = 1:10000, y = 20), rule = 2),
  DOCloss = 0,
  approx_salt = approxfun(data.frame(x = 1:10000, y = 20), rule = 2),
  ksal = c(10000, 10000, 5000),
  salsens = c(0,0,1)
)

pars_gnpz <- list(
  n.in = 5,
  ninfun1 = approxfun(data.frame(x = 1:10000, y = 0.5)),
  ninfun2 = approxfun(data.frame(x = 1:10000, y = 1)),
  n.out = 0.2,
  x.in = c(0, 0, 0, 0, 0, 0, 0, 0.02, 0, 0, 0),
  kN = rbind(c(1, 1.5, .67, 6, 2, 5, 1, 1, 1,1,1)*30.97,
             c(1, 3.5, 3.48, 3.63, 3.55, 3.62, 1, 1, 1,1,1)*14),  # allometric kN
  n.frac = rbind(rep(0.009, 11),
                 rep(0.15, 11)),
  n.frac = 0.009,
  mineralization = 0.00,
  prod = c(0,1,1,1,1,1,1,0,0,0,0),
  r = c(0, 0.6, 1.09, 0.95, 1, 1.38, 0.9, 0, 0, 0, 0),
  x.i = c(0,0,0,0,0,0,0,.012,.07,.07,.08),
  s = 0.2,
  q = 1.2,
  Ks = 1, 
  fa = 0.4,
  fm = 0.1,
  yij = adj_gnpz * 10,
  B0 = adj_gnpz * 1500,
  eij = adj_gnpz * 0.66,
  d = adj_gnpz*rbeta(length(adj_gnpz), 1,2),
  w = apply(adj_gnpz, 2, function(x){if(sum(x) > 0){x/sum(x)}else{rep(0, length(x))}}),#w2.a,
  p = get_ppar(adj_gnpz),
  m = c(0,0,0,0,0,0,0,0,0.05,0.05,0.03), #c(0,0,0,0,0,0,0,0,0.03,0.03,0.01)
  cij = diag(rep(1.8,11)),#1.8,
  h = 10,
  approx_therm = approxfun(data.frame(x = c(1, 10000), y = c(10, 10)), rule = 2),
  a.phy = 1e5,
  a.bg = 0.35, 
  I0 = c(1,25,15,20,25,25,25,1,1,1,1),
  approx_light = approxfun(data.frame(x = 1:10000, y = 150), rule = 2),
  Q10 = 3,#c(1,1,1,1,1,1,1,3,3,3,3),
  T0 = 10,
  T0.phyt = c(1,25,28,20,27,25,25,1,1,1,1),
  tsteep = c(10,150,150,150,150,150,150,10,10,10,10),
  maxt = c(40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40),
  approx_temp = approxfun(data.frame(x = 1:10000, y = 20), rule = 2),
  DOCloss = 0,
  approx_salt = approxfun(data.frame(x = 1:10000, y = 20), rule = 2),
  ksal = c(rep(10000,7), rep(5000, 4)),
  salsens = c(rep(0,7), rep(1, 4))
)


pars_gnpzf <- list(
  n.in = 5,
  ninfun1 = approxfun(data.frame(x = 1:10000, y = 0.5)),
  ninfun2 = approxfun(data.frame(x = 1:10000, y = 1)),
  n.out = 0.2,
  x.in = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  kN = rbind(c(1, 9.4, 4.37, 41.1, 13.2, 34.6, 62, 1, 1, 1,1,1),
             c(1, 25.2, 11.8, 107, 35, 90.2, 160, 1, 1, 1,1,1)),  # allometric kN
  #kN = c(1,5,3,5,4,7,5,1,1,1,1, 1),        # arbitrary kN
  n.frac = rbind(rep(0.009, 12),
                 rep(0.15, 12)),
  mineralization = 0.00,
  prod = c(0,1,1,1,1,1,1,0,0,0,0,0),
  r = c(0, 0.6, 1.09, 0.95, 1, 1.38, 0.9, 0, 0, 0, 0, 0),
  x.i = c(0,0,0,0,0,0,0,.012,.07,.07,.08,.05),
  s = 0.2,
  q = 1.2,
  Ks = 1, 
  fa = 0.4,
  fm = 0.1,
  yij = adj_gnpzf * 10,
  B0 = adj_gnpzf * 1500,
  eij = adj_gnpzf * 0.66,
  d = adj_gnpzf*rbeta(length(adj_gnpzf), 1,2),
  w = apply(adj_gnpzf, 2, function(x){if(sum(x) > 0){x/sum(x)}else{rep(0, length(x))}}),#w2.a,
  p = get_ppar(adj_gnpzf),
  m = c(0,0,0,0,0,0,0,0,0.05,0.05,0.03,0.001), #c(0,0,0,0,0,0,0,0,0.03,0.03,0.01)
  cij = diag(rep(1.8,12)),#1.8,
  h = 10,
  approx_therm = approxfun(data.frame(x = c(1, 10000), y = c(10, 10)), rule = 2),
  a.phy = 1e5,
  a.bg = 0.35, 
  I0 = c(1,25,15,20,25,25,25,1,1,1,1,1),
  approx_light = approxfun(data.frame(x = 1:10000, y = 150), rule = 2),
  Q10 = 3,#c(1,1,1,1,1,1,1,3,3,3,3),
  T0 = 10,
  T0.phyt = c(1,25,28,20,27,25,25,1,1,1,1,1),
  tsteep = c(10,150,150,150,150,150,150,10,10,10,10,10),
  maxt = c(40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40),
  approx_temp = approxfun(data.frame(x = 1:10000, y = 20), rule = 2),
  DOCloss = 0,
  approx_salt = approxfun(data.frame(x = 1:10000, y = 20), rule = 2),
  ksal = c(rep(10000,7), rep(5000, 4), 10000),
  salsens = c(rep(0,7), rep(1, 4), 0)
)

paramlst <- c(list(pars_npz, pars_gnpz, pars_gnpzf))#, plist)

#slist <- readRDS("data/smWebStates.rds")

s1 <- c(N1 = 1.5, N2 = 150, D = 100, P = 327.6, Z = 283)
s2 <- c(N1 = 1.5, N2 = 150, DOC = 100, Nanno = 50, Blue = 25, Brown = 164, Green = 20.5, Red = 36.7,
           Flag = 31.4, Rot = 3, Clad = 60, Cal = 100, Cyc = 120)
s3 <- c(N1 = 1.5, N2 = 150, DOC = 100, Nanno = 50, Blue = 25, Brown = 164, Green = 20.5, Red = 36.7,
           Flag = 31.4, Rot = 3, Clad = 60, Cal = 100, Cyc = 120, Fsh = 200)
allstates <- c(list(s1, s2, s3))#, slist)

select_foodweb <- function(x){
  return(list(mat = adj[[x]], params = paramlst[[x]], state = allstates[[x]]))
}
