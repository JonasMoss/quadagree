x <- dat.zapf2016
x <- rbind(x, x, x, x, x, x, x, x, x, x)
x <- rbind(x, x, x, x, x, x, x, x, x, x)
nrow(x) #5000

profvis::profvis(fleissci(x, boots = TRUE))
