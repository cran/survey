library(survey)


# Load dataset
data(scd)

# Create rep weights
repweights <- 2 * cbind(
  c(1,0,1,0,1,0), 
  c(1,0,0,1,0,1), 
  c(0,1,1,0,0,1),
  c(0,1,0,1,1,0)
)

scdrep <- svrepdesign(
  data = scd, 
  type = "BRR", 
  repweights = repweights, 
  combined.weights = FALSE,
  degf = 2
)


df1<-scdrep$degf

scdsub<-subset(
  x = scdrep,
  arrests >= 200
)


stopifnot(identical(scdrep$degf, scdsub$degf))
