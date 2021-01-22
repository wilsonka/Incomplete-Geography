################################################################################
## Binomial Model
################################################################################

logit <- function(p) {
  log(p/(1-p))
}
expit <- function(x) {
  exp(x)/(1+exp(x))
}

## Fit INLA model

runBinModel <- function(A, dat, X, spde, prevresult=NULL, design=NULL, ...) {
  formula <- y ~ -1 + X + f(i, model=spde)
  
  stk.dat <- inla.stack(data=dat,
                        A=list(A, 1), 
                        effects=list(i=1:spde$n.spde, 
                                     X=X),
                        tag='dat')
  
  if(!is.null(design)) {
    if(!is.null(prevresult)) {
      res <- inla(formula, family="binomial", Ntrials=N,
                  control.compute=list(dic=TRUE, mlik=T, config=T),
                  data=inla.stack.data(stk.dat),
                  control.predictor=list(A=inla.stack.A(stk.dat)#, compute=TRUE
                                         ),
                  control.fixed=list(mean=list(b0=0),
                                     prec=list(b0=1/100)),
                  control.inla=list(print.joint.hyper=T, 
                                    int.strategy="user",
                                    int.design=design,...),
                  control.mode=list(result=prevresult, restart = TRUE),
                  num.threads=1)
    } else {
      res <- inla(formula, family="binomial", Ntrials=N,
                  control.compute=list(dic=TRUE, mlik=T, config=T),
                  data=inla.stack.data(stk.dat),
                  control.predictor=list(A=inla.stack.A(stk.dat)#, compute=TRUE
                                         ),
                  control.fixed=list(mean=list(b0=0),
                                     prec=list(b0=1/100)),
                  control.inla=list(print.joint.hyper=T,
                                    int.strategy="user",
                                    int.design=design,...),
                  num.threads=1)
    }
  } else {
    res <- inla(formula, family="binomial", Ntrials=N,
                control.compute=list(dic=TRUE, mlik=T, config=T),
                data=inla.stack.data(stk.dat),
                control.predictor=list(A=inla.stack.A(stk.dat)#, compute=TRUE
                                       ),
                control.fixed=list(mean=list(b0=0),
                                   prec=list(b0=1/100)),
                control.inla=list(print.joint.hyper=T,...),
                num.threads=1)
  }
  
  return(res)
}

# Get samples from INLA model
getSamples <- function(res, nsamps) {
  samples <- inla.posterior.sample(nsamps, res)

  betas <- sapply(samples, function(x) {
    x$latent[res$misc$configs$contents$start[4]:(length(x$latent))] # assumes the last several (4+) are for the fixed effects
  })
  
  ws <- sapply(samples, function(x) {
    x$latent[res$misc$configs$contents$start[3]:(res$misc$configs$contents$start[4]-1)]
  })

  thetas <- sapply(samples, function(x) {
    x$hyperpar
  })


  list(
    beta=betas,
    w=ws,
    theta=thetas
  )
}

