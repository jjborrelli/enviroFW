library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(broom)
library(nls.multstart)

zoop <- read.csv("dataverse_files/Collated GLEON salt experiment_mesocosm data.csv")
zoop

zoop %>% 
  #filter(Lake == "Sturgeon") %>% 
  select(Mean_Cl_mg_L, Chla_ug_L, cladoceran_L) %>%
  mutate(cladoceran_L_mort = 1-cladoceran_L/max(cladoceran_L)) %>% 
  ggplot(aes(x = Mean_Cl_mg_L, y = cladoceran_L)) + geom_point()

zoop %>% 
  ggplot(aes(x = Mean_Cl_mg_L, y = Chla_ug_L)) + geom_point()



logistic <- function(x, x0.5, r){
  res <- 1 / (1 + exp(-r*(x-x0.5)))
  res[res<0] <- 0
  res
}


logistic2 <- function(x, ymax, x0.5, r){
  res <- ymax / (1 + exp(-r*(x-x0.5)))
  #res[res<0] <- 0
  res
}

unique(zoop$Lake)
testA <- zoop %>% filter(Lake == "Sturgeon") %>% 
  select(Mean_Cl_mg_L, Chla_ug_L, cladoceran_L) %>%
  filter(!is.na(cladoceran_L)) %>% 
  mutate(cladoceran_L_mort = 1 - cladoceran_L/max(cladoceran_L))

ggplot(testA, aes(x = Mean_Cl_mg_L, y = cladoceran_L_mort)) + geom_point()

fit <- nls_multstart(cladoceran_L_mort~logistic2(Mean_Cl_mg_L, ymax, x0.5, r), 
                     data = testA, iter = 1000, 
                     start_lower = c((max(testA$cladoceran_L)/2),
                                     (max(testA$cladoceran_L)/2) - (max(testA$cladoceran_L)*0.2),-.1),
                     start_upper = c((max(testA$cladoceran_L)),
                                     (max(testA$cladoceran_L)/2) + (max(testA$cladoceran_L)*0.2), 0.1))
summary(fit)
ci2 <- nlstools::confint2(fit)
ci2 


# testA %>% 
#   ggplot(aes(x = Mean_Cl_mg_L, y = cladoceran_L_mort)) + geom_point() + 
#   geom_function(fun = logistic, 
#                 args = list(x0.5 = coef(fit)[[1]],
#                             r = coef(fit)[[2]]), color = "blue", size = 1.2) + 
#   geom_function(fun = logistic, 
#                 args = list(x0.5 = ci2[1,1],
#                             r = ci2[2,2]), color = "green4",
#                 size = 1.2, lty = 2) + 
#   geom_function(fun = logistic, 
#                 args = list(x0.5 = ci2[1,2],
#                             r = ci2[2,1]), color = "green4", 
#                 size = 1.2, lty = 2)


zoop %>% #filter(Lake == "KBS Reservoir") %>% 
  select(Mean_Cl_mg_L, Acanthocyclops:Chironomid) %>%
  gather(key = taxa, value = value, Acanthocyclops:Chironomid) %>% 
  separate(taxa, into = c("gen", "spp"), sep = "_") %>% 
  select(-spp)

test <- zoop %>% #filter(Lake == "KBS Reservoir") %>% 
  select(Mean_Cl_mg_L, Acanthocyclops:Chironomid) %>%
  gather(key = taxa, value = value, Acanthocyclops:Chironomid) %>%
  separate(taxa, into = c("gen", "spp"), sep = "_") %>% 
  select(-spp) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(!is.na(value)) %>% 
  group_by(gen) %>% 
  #mutate(value = value/max(value)) %>% 
  nest() %>% 
  mutate(fit = purrr::map(data, ~ nls_multstart(value~logistic(Mean_Cl_mg_L, ymax, x0.5, r), 
                                                data = .x, iter = 1000, 
                                                start_lower = c(min(.x$value) + (max(.x$value)/2),
                                                                (max(.x$value)/2) - (max(.x$value)*0.2),-.1),
                                                start_upper = c(max(.x$value),
                                                                (max(.x$value)/2) + (max(.x$value)*0.2), 0.1),
                                                supp_errors = "Y")))

test %>% 
  mutate(aug = map(fit, augment)) %>% 
  unnest(aug) %>% 
  ggplot(aes(x = Mean_Cl_mg_L, y = value)) + geom_point() + 
  geom_line(aes(y = .fitted), color = "blue") +
  facet_wrap(~gen, scales = "free_y")


test %>% 
  mutate(par = map(fit, tidy)) %>% 
  unnest(par) %>% filter(term == "r") %>% arrange(desc(estimate)) %>% 
  print(n=150)


test %>% 
  mutate(par = map(fit, tidy)) %>% 
  unnest(par) %>% filter(term == "x0.5") %>% arrange(desc(estimate)) %>% 
  print(n=150)


params <- test %>%
  mutate(par = map(fit, tidy)) %>% 
  unnest(par)

# get confidence intervals
CI <- test %>% 
  mutate(ci = map(fit, ~ data.frame(nlstools::confint2(.x)))) %>%
  unnest(ci) %>% #colnames()
  rename(conf.low = X2.5.., conf.high = X97.5..) %>%
  group_by(gen) %>%
  mutate(term = c('ymax', 'x0.5', 'r')) %>%
  ungroup() 

params <- params %>% 
  left_join(CI, by = c("gen", "term", "data", "fit"))
params

params %>% filter(term == "x0.5", conf.high < 5000) %>% 
  ggplot(aes(x = gen, y = estimate)) + 
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) + 
  geom_point() + 
  coord_flip()




test <- zoop %>% #filter(Lake == "KBS Reservoir") %>% 
  select(Lake, Mean_Cl_mg_L, cladoceran_L:rotifer_L) %>%
  gather(key = taxa, value = value, cladoceran_L:rotifer_L) %>%
  mutate(value = as.numeric(value)) %>% 
  filter(!is.na(value)) %>% 
  group_by(taxa, Lake) %>% 
  #mutate(value = 1-value/max(value)) %>% 
  nest() %>% 
  mutate(fit = purrr::map(data, ~ nls_multstart(value~logistic2(Mean_Cl_mg_L, ymax, x0.5, r), 
                                                data = .x, iter = 1000, 
                                                start_lower = c((max(.x$value)/2),
                                                                (max(.x$value)/2) - (max(.x$value)*0.2),-.1),
                                                start_upper = c((max(.x$value)),
                                                                (max(.x$value)/2) + (max(.x$value)*0.2), 0.1),
                                                supp_errors = "Y")))

modconv <- test %>% mutate(gl = map(fit, glance)) %>% unnest(gl) %>% 
  filter(isConv)

test %>% inner_join(modconv, by = c("Lake","taxa","data", "fit")) %>% 
  filter(taxa == "cladoceran_L") %>% 
  mutate(aug = map(fit, augment)) %>% 
  unnest(aug) %>% 
  ggplot(aes(x = Mean_Cl_mg_L, y = value)) + geom_point() + 
  geom_line(aes(y = .fitted, group = Lake), color = "blue") +
  facet_wrap(~Lake, scales = "free_y") #+ scale_x_log10()


params <- test %>% inner_join(modconv, by = c("Lake","taxa","data", "fit")) %>% 
  mutate(par = map(fit, tidy)) %>% 
  unnest(par) %>% select(-colnames(modconv)[c(5:8, 10:13)])

# get confidence intervals
CI <- test %>% inner_join(modconv, by = c("Lake","taxa","data", "fit")) %>% 
  mutate(ci = map(fit, ~ data.frame(nlstools::confint2(.x)))) %>%
  unnest(ci) %>% #colnames()
  rename(conf.low = X2.5.., conf.high = X97.5..) %>%
  group_by(Lake, taxa) %>%
  mutate(term = c('ymax','x0.5', 'r')) %>%
  ungroup() %>% select(-colnames(modconv)[5:13])

params <- params %>% 
  left_join(CI, by = c("taxa", "term", "data", "fit", "Lake"))
params

params %>% filter(term == "r") %>% 
  ggplot(aes(x = taxa, y = estimate)) + 
  geom_linerange(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  geom_point() + 
  coord_flip() + facet_wrap(~Lake, scales = "free")

params %>% group_by(taxa, term) %>% summarize(value = median(estimate)) %>% spread(key = term, value = value)


params %>% filter(term == "x0.5") %>% 
  ggplot(aes(x = Lake, y = estimate)) + geom_point() + 
  coord_flip() + facet_wrap(~taxa)


fit

fit.co <- list()
for(i in 101:1000){
  newdf <- data.frame(cladoceran_L_mort = predict(fit) + sample(residuals(fit)),
                      Mean_Cl_mg_L = testA$Mean_Cl_mg_L)
  fit.new <- nls_multstart(cladoceran_L_mort~logistic(Mean_Cl_mg_L, x0.5, r), 
                           data = newdf, iter = 1000, 
                           start_lower = c((max(testA$cladoceran_L)/2) - (max(testA$cladoceran_L)*0.2),-.1),
                           start_upper = c((max(testA$cladoceran_L)/2) + (max(testA$cladoceran_L)*0.2), 0.1))
  fit.co[[i]] <-coef(fit.new)
}

apply(do.call(rbind, fit.co), 2, quantile, probs = c(0.025, 0.975)) 



testA %>% 
  ggplot(aes(x = Mean_Cl_mg_L, y = cladoceran_L_mort)) + geom_point() + 
  geom_function(fun = logistic, 
                args = list(x0.5 = coef(fit)[[1]],
                            r = coef(fit)[[2]]), color = "blue", size = 1.2) + 
  geom_function(fun = logistic, 
                args = list(x0.5 = -337,
                            r = 0.007675), color = "green4",
                size = 1.2, lty = 2) + 
  geom_function(fun = logistic, 
                args = list(x0.5 = 125,
                            r = 0.001474765), color = "green4", 
                size = 1.2, lty = 2)