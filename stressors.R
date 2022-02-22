

#### TEMPERATURE ####
tmpdf_const <- data.frame(x = 1:10000, y = 20)
tmpdf_1 <- data.frame(x = 1:10000, y = 1*sin(((1:10000) - 0)*(.017))+20)
tmpdf_5 <- data.frame(x = 1:10000, y = 5*sin(((1:10000) - 0)*(.017))+20)

tmpdf_1_stoch <- data.frame(x = 1:10000, y = abs(rnorm(10000, 1*sin(((1:10000) - 0)*(.017))+20, 1)))
plot(tmpdf_1_stoch[1:3650,])

tmpdf_5_stoch <- data.frame(x = 1:10000, y = abs(rnorm(10000, 5*sin(((1:10000) - 0)*(.017))+20, 1)))
plot(tmpdf_5_stoch[1:3650,])

#### LIGHT ####
lmu <- sample(200:300, 10000, replace = TRUE)
lmu2 <- sample(25:50, 10000, replace = TRUE)
lgtdf <- data.frame(x = 1:10000, y = abs(rnorm(10000, lmu2*sin(((1:10000) - 0)*(.017))+lmu, 10)))

lgtdf_const_low <- data.frame(x = 1:10000, y = 150)
lgtdf_const_hi <- data.frame(x = 1:10000, y = 250)
yr1 <- (1:720 + ymd("2020-01-01"))
lgtdf_40 <- data.frame(x = 1:10000, y = 40*sin(((1:10000) - 0)*(.013))+200)
plot(lgtdf_40[1:720,2] ~ yr1)
lgtdf_60 <- data.frame(x = 1:10000, y = 60*sin(((1:10000) - 0)*(.013))+200)
plot(lgtdf_60[1:3650,])

lgtdf_40_stoch <- data.frame(x = 1:10000, 
                            y = abs(rnorm(10000, 40*sin(((1:10000) - 0)*(.017))+200, 10)))
plot(lgtdf_40_stoch[1:3650,])
lgtdf_60_stoch <- data.frame(x = 1:10000, 
                            y = abs(rnorm(10000, 60*sin(((1:10000) - 0)*(.017))+200, 10)))
plot(lgtdf_60_stoch[1:3650,])

#### NUTRIENTS ####
nmu <- (runif(10000, 0.5, 2))
ndf <- data.frame(x = 1:10000, 
                  y = runif(10000,0,1)*abs(rnorm(10000, sin(((1:10000) - 0)*(.017))+nmu, 0.25)))

plot(ndf[1:3650,])

ndf_const_low <- data.frame(x = 1:10000, y = 0.5)
ndf_const_hi <- data.frame(x = 1:10000, y = 2.5)


ndf_1 <- data.frame(x = 1:10000, y = 1+(1*sin(((1:10000) - 0)*(.017))))
plot(ndf_1[1:365,])
ndf_2 <- data.frame(x = 1:10000, y = abs(2*sin(((1:10000) - 0)*(.017))+sin(((1:10000) - 0)*(.017*2))))
plot(ndf_2[1:365,])
