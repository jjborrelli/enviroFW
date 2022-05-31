

tmpdf <- data.frame(x = 1:10000, y = 20)
lgtdf <- data.frame(x = 1:10000, y = 150)
ndf_1 <- data.frame(x = 1:10000, y = 1)
ndf_2 <- data.frame(x = 1:10000, y = 5)
saldf <- data.frame(x = 1:10000, y = 15)

pars_npz$ninfun1 <- approxfun(ndf_1)
pars_npz$ninfun2 <- approxfun(ndf_2)
pars_npz$approx_temp <- approxfun(tmpdf)
pars_npz$approx_light <- approxfun(lgtdf)
pars_npz$approx_salt <- approxfun(saldf)

# run dynamics in constant environment

# pars_npz$ksal
pars_npz$salsens <- c(0,0,1)
pars_npz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 0))
out <- ode(s1, 1:1000, atn_2N, pars_npz, events = list(func = goExtinct, time = 1:3000))

pars_npz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 50))
out1 <- ode(s1, 1:1000, atn_2N, pars_npz, events = list(func = goExtinct, time = 1:3000))

pars_npz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 100))
out2 <- ode(s1, 1:1000, atn_2N, pars_npz, events = list(func = goExtinct, time = 1:3000))

pars_npz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 150))
out3 <- ode(s1, 1:1000, atn_2N, pars_npz, events = list(func = goExtinct, time = 1:3000))

o1 <- reshape2::melt(out) %>% 
  mutate(sal = 0)
o2 <- reshape2::melt(out1) %>% 
  mutate(sal = 50)
o3 <- reshape2::melt(out2) %>% 
  mutate(sal = 100)
o4 <- reshape2::melt(out3) %>% 
  mutate(sal = 150)

do.call(rbind, list(o1, o2, o3, o4)) %>% 
  filter(Var2 != "time") %>% 
  #filter(Var1 %in% c(750:1000)) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(sal))) + geom_line() + 
  facet_wrap(~Var2, scales = "free_y") + scale_y_log10()


# run dynamics in fluctuating  environment

pars_npz$ninfun1 <- approxfun(data.frame(x = 1:10000, y = 1+(1*sin(((1:10000) - 0)*(.017)))))
pars_npz$ninfun2 <- approxfun(data.frame(x = 1:10000, y = abs(2*sin(((1:10000) - 0)*(.017))+sin(((1:10000) - 0)*(.017*2)))))
pars_npz$approx_temp <- approxfun(data.frame(x = 1:10000, y = 5*sin(((1:10000) - 0)*(.017))+20))
# pars_npz$ksal
pars_npz$salsens <- c(0,0,1.2)
pars_npz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 0))
out <- ode(s1, 1:1000, atn_2N, pars_npz, events = list(func = goExtinct, time = 1:3000))

pars_npz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 50))
out1 <- ode(s1, 1:1000, atn_2N, pars_npz, events = list(func = goExtinct, time = 1:3000))

pars_npz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 100))
out2 <- ode(s1, 1:1000, atn_2N, pars_npz, events = list(func = goExtinct, time = 1:3000))

pars_npz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 150))
out3 <- ode(s1, 1:1000, atn_2N, pars_npz, events = list(func = goExtinct, time = 1:3000))

o1 <- reshape2::melt(out) %>% 
  mutate(sal = 0)
o2 <- reshape2::melt(out1) %>% 
  mutate(sal = 50)
o3 <- reshape2::melt(out2) %>% 
  mutate(sal = 100)
o4 <- reshape2::melt(out3) %>% 
  mutate(sal = 150)

p1 <- do.call(rbind, list(o1, o2, o3, o4)) %>% 
  filter(Var2 != "time") %>% 
  #filter(Var1 %in% c(750:1000)) %>%
  #filter(Var2 %in% c("P", "Z")) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(sal))) + geom_line() + 
  facet_wrap(~Var2, scales = "free_y", nrow = 5)
p2 <- ggplot(data.frame(x = 1:1000, y = 5*sin(((1:1000) - 0)*(.017))+20)) + 
  geom_line(aes(x = x , y = y))

p1/p2 + plot_layout(ncol = 1, heights = c(5,1))



## GNPZ

tmpdf <- data.frame(x = 1:10000, y = 20)
lgtdf <- data.frame(x = 1:10000, y = 150)
ndf_1 <- data.frame(x = 1:10000, y = 3)
ndf_2 <- data.frame(x = 1:10000, y = 6)
saldf <- data.frame(x = 1:10000, y = 15)

pars_gnpz$ninfun1 <- approxfun(ndf_1)
pars_gnpz$ninfun2 <- approxfun(ndf_2)
pars_gnpz$approx_temp <- approxfun(tmpdf)
pars_gnpz$approx_light <- approxfun(lgtdf)
pars_gnpz$approx_salt <- approxfun(saldf)

pars_gnpz$salsens <- c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0)

pars_gnpz$approx_temp <- approxfun(data.frame(x = 1:10000, y = 5*sin(((1:10000) - 0)*(.017))+20))
pars_gnpz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 0))
out <- ode(s2, 1:1000, atn_2N, pars_gnpz, events = list(func = goExtinct, time = 1:1000))
tail(out)

pars_gnpz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 50))
out1 <- ode(s2, 1:1000, atn_2N, pars_gnpz, events = list(func = goExtinct, time = 1:1000))
tail(out1)

pars_gnpz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 100))
out2 <- ode(s2, 1:1000, atn_2N, pars_gnpz, events = list(func = goExtinct, time = 1:1000))
tail(out2)

pars_gnpz$approx_salt <- approxfun(data.frame(x = 1:10000, y = 150))
out3 <- ode(s2, 1:1000, atn_2N, pars_gnpz, events = list(func = goExtinct, time = 1:1000))
tail(out3)

o1 <- reshape2::melt(out) %>% 
  mutate(sal = 0)
o2 <- reshape2::melt(out1) %>% 
  mutate(sal = 50)
o3 <- reshape2::melt(out2) %>% 
  mutate(sal = 100)
o4 <- reshape2::melt(out3) %>% 
  mutate(sal = 150)

p1 <- do.call(rbind, list(o1, o2, o3, o4)) %>% 
  filter(Var2 != "time") %>% 
  #filter(Var1 %in% c(750:1000)) %>%
  #filter(Var2 %in% c("P", "Z")) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(sal))) + geom_line() + 
  facet_wrap(~Var2, scales = "free_y")
p1




## GNPZF

tmpdf <- data.frame(x = 1:10000, y = 20)
lgtdf <- data.frame(x = 1:10000, y = 150)
ndf_1 <- data.frame(x = 1:10000, y = 30)
ndf_2 <- data.frame(x = 1:10000, y = 60)
saldf <- data.frame(x = 1:10000, y = 15)

pars_gnpzf$ninfun1 <- approxfun(ndf_1)
pars_gnpzf$ninfun2 <- approxfun(ndf_2)
pars_gnpzf$approx_temp <- approxfun(tmpdf)
pars_gnpzf$approx_light <- approxfun(lgtdf)
pars_gnpzf$approx_salt <- approxfun(saldf)

pars_gnpzf$salsens <- c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0)

#pars_gnpzf$approx_temp <- approxfun(data.frame(x = 1:10000, y = 5*sin(((1:10000) - 0)*(.017))+20))
pars_gnpzf$approx_salt <- approxfun(data.frame(x = 1:10000, y = 0))
out <- ode(s3, 1:1000, atn_2N, pars_gnpzf, events = list(func = goExtinct, time = 1:1000))
tail(out)

pars_gnpzf$approx_salt <- approxfun(data.frame(x = 1:10000, y = 50))
out1 <- ode(s3, 1:1000, atn_2N, pars_gnpzf, events = list(func = goExtinct, time = 1:1000))
tail(out1)

pars_gnpzf$approx_salt <- approxfun(data.frame(x = 1:10000, y = 100))
out2 <- ode(s3, 1:1000, atn_2N, pars_gnpzf, events = list(func = goExtinct, time = 1:1000))
tail(out2)

pars_gnpzf$approx_salt <- approxfun(data.frame(x = 1:10000, y = 150))
out3 <- ode(s3, 1:1000, atn_2N, pars_gnpzf, events = list(func = goExtinct, time = 1:1000))
tail(out3)

pars_gnpzf$approx_salt <- approxfun(data.frame(x = 1:10000, y = 550))
out4 <- ode(s3, 1:1000, atn_2N, pars_gnpzf, events = list(func = goExtinct, time = 1:1000))
tail(out4)


o1 <- reshape2::melt(out) %>% 
  mutate(sal = 0)
o2 <- reshape2::melt(out1) %>% 
  mutate(sal = 50)
o3 <- reshape2::melt(out2) %>% 
  mutate(sal = 100)
o4 <- reshape2::melt(out3) %>% 
  mutate(sal = 150)
o5 <- reshape2::melt(out4) %>% 
  mutate(sal = 550)

p1 <- do.call(rbind, list(o1, o2, o3, o4, o5)) %>% 
  filter(Var2 != "time") %>% 
  #filter(Var1 %in% c(750:1000)) %>%
  #filter(Var2 %in% c("P", "Z")) %>% 
  ggplot(aes(x = Var1, y = value, color = factor(sal))) + geom_line(size = 1.2) + 
  facet_wrap(~Var2, scales = "free_y")
p1
