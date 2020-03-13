## Load in the sintering data
alumina <- read.csv("inst/extdata/sintering_data.csv")
alumina$Sample <- as.factor(alumina$Sample)

save(alumina, file="data/alumina.Rda")
