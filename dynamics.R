library(dplyr)
library(lubridate)
library(deSolve)
library(ggplot2)


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

temp_lim <- function(temp, Topt = 25, Tmax = 40, Tsteepness = 150){
  temp_opt <- exp(-(((Topt - temp)^2)/Tsteepness))
  penalty <- 1 - (1/(1 + exp(-temp + Tmax)))
  
  return(temp_opt * penalty)
}
Fij <- function(B0, q, d, p, B, w){
  numer <-  w * B^q
  halfsat <- B0^q
  interfere <- 0#sapply(1:length(B), function(k){sum(d[k,] * p[,k] * B[k] * B0[k,]^q)}) # d * p * B * B0^q
  res <- colSums(w * B^q)
  
  fres <- numer/t(apply(halfsat, 1, function(x) x + interfere + res))
  fres[is.nan(fres)] <- 0
  return(fres)
}
goExtinct <- function(times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0
    return(c(states))
  })
}


plist <- readRDS("data/smWebPARAMs.rds")
mlist <- readRDS("data/smWebMats.rds")
slist <- readRDS("data/smWebStates.rds")

cbind(sapply(mlist, nrow), sapply(mlist, sum))

tmu <- sample(17:22, 10000, replace = TRUE)
tmu2 <- sample(2:7, 10000, replace = TRUE)
tmpdf <- data.frame(x = 1:10000, y = abs(rnorm(10000, tmu2*sin(((1:10000) - 0)*(.017))+tmu, 1)))
lmu <- sample(200:300, 10000, replace = TRUE)
lmu2 <- sample(25:50, 10000, replace = TRUE)
lgtdf <- data.frame(x = 1:10000, y = abs(rnorm(10000, lmu2*sin(((1:10000) - 0)*(.017))+lmu, 10)))
nmu <- (runif(10000, 0.5, 2))
ndf <- data.frame(x = 1:10000, 
                  y = runif(10000,0,1)*abs(rnorm(10000, sin(((1:10000) - 0)*(.017))+nmu, 0.25)))

t1 <- Sys.time()
dlist <- lapply(1:6, function(x){
  plist[[x]]$ninfun <- approxfun(ndf)
  plist[[x]]$approx_therm <- approxfun(data.frame(x = 1:10000, y = rep(7, 10000)))
  plist[[x]]$approx_light <- approxfun(lgtdf)
  plist[[x]]$approx_temp <- approxfun(tmpdf)
  
  out <- ode(slist[[x]], times = 1:3650, func = atn, parms = plist[[x]], 
             events = list(func = goExtinct, time = 1:3650))
  odf <- reshape2::melt(out[,2:ncol(out)])
  odf$config <- x
  return(odf)
})
t2 <- Sys.time()
t2 - t1

dlist <- do.call(rbind, dlist)

dlist %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!Var2 %in% c("L", "tem")) %>% 
  mutate(grp = case_when(Var2 == "blue" ~ "phyto", Var2 == "brown" ~ "phyto",
                         Var2 == "green" ~ "phyto", Var2 == "red" ~ "phyto",
                         Var2 == "flag" ~ "prot", Var2 == "rot" ~ "zoop",
                         Var2 == "clad" ~ "zoop", Var2 == "cyc" ~ "zoop",
                         Var2 == "cal" ~ "zoop", Var2 == "benthFilt" ~ "macroinv",
                         Var2 == "benthDep" ~ "macroinv", Var2 == "benthPred" ~ "macroinv",
                         Var2 == "smFish" ~ "fish", Var2 == "lrgFish" ~ "fish", TRUE ~ Var2)) %>% 
  group_by(Var1, grp, config) %>% summarize(value = sum(value)) %>% 
  mutate(grp = factor(grp, levels = c("L", "tem", "N", "detritus", "plant", "phyto", "prot", "zoop",
                                      "macroinv", "fish"))) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(config))) + geom_line(size = 1.2) + 
  #scale_y_log10() +
  scale_color_viridis_d(end = 0.8) + 
  facet_wrap(~grp, scales = "free_y", nrow = 2) + 
  theme_bw() + 
  labs(color = "")


dlist %>% 
  filter(Var2 %in% c("blue", "brown", "green", "red", "rot", "clad", "cal", "cyc")) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(config), group = config)) +
  geom_line(size = 1.2) +
  facet_wrap(~Var2, scales = "free_y", nrow = 2) + 
  scale_color_viridis_d(end = 0.8) + scale_y_log10() + 
  labs(color = "")



dlist %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!Var2 %in% c("L", "tem")) %>% 
  mutate(grp = case_when(Var2 == "blue" ~ "phyto", Var2 == "brown" ~ "phyto",
                         Var2 == "green" ~ "phyto", Var2 == "red" ~ "phyto",
                         Var2 == "flag" ~ "prot", Var2 == "rot" ~ "zoop",
                         Var2 == "clad" ~ "zoop", Var2 == "cyc" ~ "zoop",
                         Var2 == "cal" ~ "zoop", Var2 == "benthFilt" ~ "macroinv",
                         Var2 == "benthDep" ~ "macroinv", Var2 == "benthPred" ~ "macroinv",
                         Var2 == "smFish" ~ "fish", Var2 == "lrgFish" ~ "fish", TRUE ~ Var2)) %>% 
  group_by(Var1, grp, config) %>% summarize(value = sum(value)) %>% 
  mutate(grp = factor(grp, levels = c("L", "tem", "N", "detritus", "plant", "phyto", "prot", "zoop",
                                      "macroinv", "fish"))) %>% 
  filter(grp %in% c("phyto","zoop")) %>% 
  ggplot(aes(x = Var1+ymd("2015-01-01"), y = value, color = grp)) + geom_line(size = 1.2) + 
  scale_y_log10() +
  scale_color_viridis_d(end = 0.8) + 
  facet_wrap(~config, scales = "free_y", nrow = 2) + 
  theme_bw() + 
  labs(color = "", x = "Date", y = "Biomass")


##################################################
# Vary temp
tmu <- sample(17:22, 10000, replace = TRUE)
tmu2 <- sample(2:7, 10000, replace = TRUE)
tmpdf <- data.frame(x = 1:10000, y = abs(rnorm(10000, tmu2*sin(((1:10000) - 0)*(.017))+tmu, 1)))
lgtdf <- data.frame(x = 1:10000, y = 250)
ndf <- data.frame(x = 1:10000, y = 1)

t1 <- Sys.time()
dlistTEMP <- lapply(1:6, function(x){
  plist[[x]]$ninfun <- approxfun(ndf)
  plist[[x]]$approx_therm <- approxfun(data.frame(x = 1:10000, y = rep(7, 10000)))
  plist[[x]]$approx_light <- approxfun(lgtdf)
  plist[[x]]$approx_temp <- approxfun(tmpdf)
  
  out <- ode(slist[[x]], times = 1:3650, func = atn, parms = plist[[x]], 
             events = list(func = goExtinct, time = 1:3650))
  odf <- reshape2::melt(out[,2:ncol(out)])
  odf$config <- x
  return(odf)
})
t2 <- Sys.time()
t2 - t1

dlistTEMP <- do.call(rbind, dlistTEMP)

dlistTEMP %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!Var2 %in% c("L", "tem")) %>% 
  mutate(grp = case_when(Var2 == "blue" ~ "phyto", Var2 == "brown" ~ "phyto",
                         Var2 == "green" ~ "phyto", Var2 == "red" ~ "phyto",
                         Var2 == "flag" ~ "prot", Var2 == "rot" ~ "zoop",
                         Var2 == "clad" ~ "zoop", Var2 == "cyc" ~ "zoop",
                         Var2 == "cal" ~ "zoop", Var2 == "benthFilt" ~ "macroinv",
                         Var2 == "benthDep" ~ "macroinv", Var2 == "benthPred" ~ "macroinv",
                         Var2 == "smFish" ~ "fish", Var2 == "lrgFish" ~ "fish", TRUE ~ Var2)) %>% 
  group_by(Var1, grp, config) %>% summarize(value = sum(value)) %>% 
  mutate(grp = factor(grp, levels = c("L", "tem", "N", "detritus", "plant", "phyto", "prot", "zoop",
                                      "macroinv", "fish"))) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(config))) + geom_line(size = 1.2) + 
  #scale_y_log10() +
  scale_color_viridis_d(end = 0.8) + 
  facet_wrap(~grp, scales = "free_y", nrow = 2) + 
  theme_bw() + 
  labs(color = "")


##################################################
# Vary light

tmpdf <- data.frame(x = 1:10000, y = 20)
lmu <- sample(200:300, 10000, replace = TRUE)
lmu2 <- sample(25:50, 10000, replace = TRUE)
lgtdf <- data.frame(x = 1:10000, y = abs(rnorm(10000, lmu2*sin(((1:10000) - 0)*(.017))+lmu, 10)))
nmu <- (runif(10000, 0.5, 2))
ndf <- data.frame(x = 1:10000, y = 1)


t1 <- Sys.time()
dlistLGHT <- lapply(1:6, function(x){
  plist[[x]]$ninfun <- approxfun(ndf)
  plist[[x]]$approx_therm <- approxfun(data.frame(x = 1:10000, y = rep(7, 10000)))
  plist[[x]]$approx_light <- approxfun(lgtdf)
  plist[[x]]$approx_temp <- approxfun(tmpdf)
  
  out <- ode(slist[[x]], times = 1:3650, func = atn, parms = plist[[x]], 
             events = list(func = goExtinct, time = 1:3650))
  odf <- reshape2::melt(out[,2:ncol(out)])
  odf$config <- x
  return(odf)
})
t2 <- Sys.time()
t2 - t1

dlistLGHT <- do.call(rbind, dlistLGHT)

dlistLGHT %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!Var2 %in% c("L", "tem")) %>% 
  mutate(grp = case_when(Var2 == "blue" ~ "phyto", Var2 == "brown" ~ "phyto",
                         Var2 == "green" ~ "phyto", Var2 == "red" ~ "phyto",
                         Var2 == "flag" ~ "prot", Var2 == "rot" ~ "zoop",
                         Var2 == "clad" ~ "zoop", Var2 == "cyc" ~ "zoop",
                         Var2 == "cal" ~ "zoop", Var2 == "benthFilt" ~ "macroinv",
                         Var2 == "benthDep" ~ "macroinv", Var2 == "benthPred" ~ "macroinv",
                         Var2 == "smFish" ~ "fish", Var2 == "lrgFish" ~ "fish", TRUE ~ Var2)) %>% 
  group_by(Var1, grp, config) %>% summarize(value = sum(value)) %>% 
  mutate(grp = factor(grp, levels = c("L", "tem", "N", "detritus", "plant", "phyto", "prot", "zoop",
                                      "macroinv", "fish"))) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(config))) + geom_line(size = 1.2) + 
  #scale_y_log10() +
  scale_color_viridis_d(end = 0.8) + 
  facet_wrap(~grp, scales = "free_y", nrow = 2) + 
  theme_bw() + 
  labs(color = "")


##################################################
# Vary N

tmpdf <- data.frame(x = 1:10000, y = 20)
lgtdf <- data.frame(x = 1:10000, y = 250)
nmu <- (runif(10000, 0.5, 2))
ndf <- data.frame(x = 1:10000, 
                  y = runif(10000,0,1)*abs(rnorm(10000, sin(((1:10000) - 0)*(.017))+nmu, 0.25)))

t1 <- Sys.time()
dlistN <- lapply(1:6, function(x){
  plist[[x]]$ninfun <- approxfun(ndf)
  plist[[x]]$approx_therm <- approxfun(data.frame(x = 1:10000, y = rep(7, 10000)))
  plist[[x]]$approx_light <- approxfun(lgtdf)
  plist[[x]]$approx_temp <- approxfun(tmpdf)
  
  out <- ode(slist[[x]], times = 1:3650, func = atn, parms = plist[[x]], 
             events = list(func = goExtinct, time = 1:3650))
  odf <- reshape2::melt(out[,2:ncol(out)])
  odf$config <- x
  return(odf)
})
t2 <- Sys.time()
t2 - t1

dlistN <- do.call(rbind, dlistN)

dlistN %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!Var2 %in% c("L", "tem")) %>% 
  mutate(grp = case_when(Var2 == "blue" ~ "phyto", Var2 == "brown" ~ "phyto",
                         Var2 == "green" ~ "phyto", Var2 == "red" ~ "phyto",
                         Var2 == "flag" ~ "prot", Var2 == "rot" ~ "zoop",
                         Var2 == "clad" ~ "zoop", Var2 == "cyc" ~ "zoop",
                         Var2 == "cal" ~ "zoop", Var2 == "benthFilt" ~ "macroinv",
                         Var2 == "benthDep" ~ "macroinv", Var2 == "benthPred" ~ "macroinv",
                         Var2 == "smFish" ~ "fish", Var2 == "lrgFish" ~ "fish", TRUE ~ Var2)) %>% 
  group_by(Var1, grp, config) %>% summarize(value = sum(value)) %>% 
  mutate(grp = factor(grp, levels = c("L", "tem", "N", "detritus", "plant", "phyto", "prot", "zoop",
                                      "macroinv", "fish"))) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(config))) + geom_line(size = 1.2) + 
  #scale_y_log10() +
  scale_color_viridis_d(end = 0.8) + 
  facet_wrap(~grp, scales = "free_y", nrow = 2) + 
  theme_bw() + 
  labs(color = "")


##################################################
# Vary N and Temp

tmu <- sample(17:22, 10000, replace = TRUE)
tmu2 <- sample(2:7, 10000, replace = TRUE)
tmpdf <- data.frame(x = 1:10000, 
                    #y = abs(rnorm(10000, 5*sin(((1:10000) - 0)*(.017))+20, 1))
                    y = 2*sin(((1:10000) - 0)*(.017))+18)
lgtdf <- data.frame(x = 1:10000, y = 250)
nmu <- (runif(10000, 0.5, 2))
ndf <- data.frame(x = 1:10000, 
                  #y = runif(10000,0,1)*abs(rnorm(10000, sin(((1:10000) - 0)*(.017))+1, 0.25))
                  y = sin(((1:10000) + 30)*(.017))+1)

t1 <- Sys.time()
dlistNT <- lapply(1:6, function(x){
  plist[[x]]$ninfun <- approxfun(ndf)
  plist[[x]]$approx_therm <- approxfun(data.frame(x = 1:10000, y = rep(7, 10000)))
  plist[[x]]$approx_light <- approxfun(lgtdf)
  plist[[x]]$approx_temp <- approxfun(tmpdf)
  
  out <- ode(slist[[x]], times = 1:3650, func = atn, parms = plist[[x]], 
             events = list(func = goExtinct, time = 1:3650))
  odf <- reshape2::melt(out[,2:ncol(out)])
  odf$config <- x
  return(odf)
})
t2 <- Sys.time()
t2 - t1

dlistNT <- do.call(rbind, dlistNT)

dlistNT %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!Var2 %in% c("L", "tem")) %>% 
  mutate(grp = case_when(Var2 == "blue" ~ "phyto", Var2 == "brown" ~ "phyto",
                         Var2 == "green" ~ "phyto", Var2 == "red" ~ "phyto",
                         Var2 == "flag" ~ "prot", Var2 == "rot" ~ "zoop",
                         Var2 == "clad" ~ "zoop", Var2 == "cyc" ~ "zoop",
                         Var2 == "cal" ~ "zoop", Var2 == "benthFilt" ~ "macroinv",
                         Var2 == "benthDep" ~ "macroinv", Var2 == "benthPred" ~ "macroinv",
                         Var2 == "smFish" ~ "fish", Var2 == "lrgFish" ~ "fish", TRUE ~ Var2)) %>% 
  group_by(Var1, grp, config) %>% summarize(value = sum(value)) %>% 
  mutate(grp = factor(grp, levels = c("L", "tem", "N", "detritus", "plant", "phyto", "prot", "zoop",
                                      "macroinv", "fish"))) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(config))) + geom_line(size = 1.2) + 
  scale_y_log10() +
  scale_color_viridis_d(end = 0.8) + 
  facet_wrap(~grp, scales = "free_y", nrow = 2) + 
  theme_bw() + 
  labs(color = "")


dlistNT %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!Var2 %in% c("L", "tem")) %>% 
  mutate(grp = case_when(Var2 == "blue" ~ "phyto", Var2 == "brown" ~ "phyto",
                         Var2 == "green" ~ "phyto", Var2 == "red" ~ "phyto",
                         Var2 == "flag" ~ "prot", Var2 == "rot" ~ "zoop",
                         Var2 == "clad" ~ "zoop", Var2 == "cyc" ~ "zoop",
                         Var2 == "cal" ~ "zoop", Var2 == "benthFilt" ~ "macroinv",
                         Var2 == "benthDep" ~ "macroinv", Var2 == "benthPred" ~ "macroinv",
                         Var2 == "smFish" ~ "fish", Var2 == "lrgFish" ~ "fish", TRUE ~ Var2)) %>% 
  filter(grp == "phyto", config %in% 3:6) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(config))) + geom_line(size = 1.2) + 
  #scale_y_log10() +
  scale_color_viridis_d(end = 0.8) + 
  facet_wrap(~Var2, scales = "free_y", nrow = 2) + 
  theme_bw() + 
  labs(color = "")



dlistNT %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!Var2 %in% c("L", "tem")) %>% 
  mutate(grp = case_when(Var2 == "blue" ~ "phyto", Var2 == "brown" ~ "phyto",
                         Var2 == "green" ~ "phyto", Var2 == "red" ~ "phyto",
                         Var2 == "flag" ~ "prot", Var2 == "rot" ~ "zoop",
                         Var2 == "clad" ~ "zoop", Var2 == "cyc" ~ "zoop",
                         Var2 == "cal" ~ "zoop", Var2 == "benthFilt" ~ "macroinv",
                         Var2 == "benthDep" ~ "macroinv", Var2 == "benthPred" ~ "macroinv",
                         Var2 == "smFish" ~ "fish", Var2 == "lrgFish" ~ "fish", TRUE ~ Var2)) %>% 
  filter(grp == "phyto", config %in% 3:6, Var2 == "red") %>% 
  mutate(grprun = paste(Var2, config, sep = "_")) %>% 
  select(Var1, grprun, value) %>% 
  tidyr::spread(key = grprun, value = value) %>% select(-Var1) %>% 
  GGally::ggpairs()



dlistNT %>% 
  mutate(Var2 = as.character(Var2)) %>%
  group_by(Var2, config) %>% 
  summarize(cv = sd(value)/mean(value)) %>% 
  tidyr::spread(key = config, value = cv) %>% 
  print(n = 24)

