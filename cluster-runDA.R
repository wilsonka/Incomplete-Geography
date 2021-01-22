args = commandArgs(trailingOnly = TRUE)
print(args)
i = eval(parse(text=args[1]))
spcov = eval(parse(text=args[2]))
jitter = eval(parse(text=args[3]))


print(i)
print(spcov)
print(jitter)
set.seed(i*1001)



library(INLA)

if(jitter) {
  if(spcov) {
    load("Data/forcluster-fitDA-spcov.rda")
  } else {
    load("Data/forcluster-fitDA.rda")
  }
  prior.loc <- prior_disp
  beta.init <- res.naive$summary.fixed$mean
  w.init <- res.naive$summary.random$i$mean
  idx.masked <- NULL
  mask <- F
} else {
  if(spcov) {
    load("Data/forcluster-fitDA-spcov-m.rda")
  } else {
    load("Data/forcluster-fitDA-m.rda")
  }
  prior.loc <- prior_mask
  beta.init <- res50$summary.fixed$mean
  w.init <- res50$summary.random$i$mean
  mask <- T
}


source("binomial-model.R")
source("run-DA.R")
output <- runDA(clusters08, masterframe, 
                prior.loc=prior.loc, X, beta.init=beta.init, 
                w.init=w.init, 
                niter=1000, nburnin=100, mesh, spde, spcov = spcov,
                mask=mask, idx.masked=idx.masked)
save(output, file=paste0("Results/runDA-jitter",jitter,"-cov",spcov,"-",i,".rda"))
