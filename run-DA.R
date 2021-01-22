################################################################################
## Complete Algorithm
################################################################################

expit <- function(x) {
  exp(x)/(1 + exp(x))
}

logit <- function(p) {
  log(p/(1-p))
}

# Gibbs sampling for location
updateloc.c <- function(y.c, N.c, pix.c, p.c, idx.c) {
  # For cluster c, do Gibbs sampling
  # y.c: number of outcome
  # n.c: number of trials
  # pix.c: length of number of possible enumeration areas
  #       prior on location
  # p.c: length of number of possible enumeration areas
  #       prob of outcome given that EA
  # idx.c: index of which EAs
  post <- p.c^y.c * (1-p.c)^(N.c-y.c) * pix.c
  post.sum <- sum(post)
  sample(idx.c, 1, prob=as.vector(post/post.sum))
}

determinep.c <- function(beta, X.c, A.c, w, Z.c=NULL) {
  # Determine p.c (the prob of outcome given that EA for cluster c)
  # beta: fixed effects (first p correspond to X.c)
  # X.c: design matrix for cluster level effects
  # A.c mapping from EAs to mesh points
  # w: spatial random effect
  # Z.c: spatial covariates (correspond to value at EA)
  
  if(!is.null(Z.c)) {
    beta.clust <- beta[1:length(X.c)]
    beta.sp <- beta[(length(X.c) + 1): length(beta)]
    expit(matrix(X.c %*% beta.clust, nrow=nrow(A.c), ncol=length(X.c)) + Z.c %*% matrix(beta.sp, ncol=1) + A.c %*% w)
  } else {
    expit(matrix(X.c %*% beta, nrow=nrow(A.c), ncol=length(X.c)) + A.c %*% w)
  }
  
}

runone <- function(clusters, masterframe, prior.loc, beta, w, X, Alist, spde, mesh, spcov=F, mask=F, idx.masked=NULL) {
  ## Update location
  idx <- sapply(1:nrow(clusters[clusters$mask,]), function(i) {
    if(mask) {
      c <- idx.masked[i]
    } else {
      c <- i
    }
    if(length(prior.loc$idxlist[[i]])==1) {
      prior.loc$idxlist[[i]]
    } else {
      if(spcov) { # assumes the masterframe contains the correct transformation 
        # and there is only one spatial covariate
        p.c <- determinep.c(beta, X[c, ], Alist[[i]], w, 
                            Z.c=masterframe$spcov[prior.loc$idxlist[[i]]])
        updateloc.c(clusters$y[c], clusters$N[c], prior.loc$prior[[i]], 
                    p.c, prior.loc$idxlist[[i]])
      } else {
        p.c <- determinep.c(beta, X[c, ], Alist[[i]], w)
        updateloc.c(clusters$y[c], clusters$N[c], prior.loc$prior[[i]], 
                    p.c, prior.loc$idxlist[[i]])
      }
    }
  })
  
  ## Update other parameters
  if(spcov) {
    formula <- ~ -1 + b0 + light
    if(mask) {
      df.inla <- data.frame(b0=rep(1, nrow(clusters)),
                            light=rep(1, nrow(clusters)))
      df.inla$light[-idx.masked] <- clusters08$spcov_correct[-idx.masked]
      df.inla$light[idx.masked] <- masterframe$spcov[idx]
      X.inla <- model.matrix(formula, df.inla)
    } else {
      X.inla <- model.matrix(formula, data.frame(b0=rep(1, nrow(clusters)),
                                                 light=masterframe$spcov[idx]))
    }
    
  } else {
    formula <- ~ -1 + b0
    X.inla <- model.matrix(formula, data.frame(b0=rep(1, nrow(clusters))))
  }
  # Step 1: Get things ready for INLA model
  if(mask) {
    locmat <- matrix(0, nrow=nrow(clusters), ncol=2)
    locmat[-idx.masked, ] <- as.matrix(clusters[-idx.masked, c("east.correct", "north.correct")])
    locmat[idx.masked, ] <- as.matrix(masterframe[idx,c("east", "north")])
  } else {
    locmat <- as.matrix(masterframe[idx,c("east", "north")])
  }
  A <- inla.spde.make.A(mesh, loc=locmat)
  # Step 2: Run INLA model
  res.inla <- runBinModel(A=A, dat=list(y=clusters$y, N=clusters$N), 
                          X=X.inla, spde=spde)
  # Step 3: Sample from INLA model
  samp <- getSamples(res.inla, nsamps=1)
  
  list(theta=samp$theta,
       beta=samp$beta,
       w=samp$w,
       loc=idx)
}

runDA <- function(clusters, masterframe, prior.loc, X, beta.init, 
                  w.init, niter, nburnin, mesh, spde, spcov=F, mask=F, idx.masked=NULL) {
  betas <- matrix(0, nrow=niter, ncol=length(beta.init))
  ws <- matrix(0, nrow=niter, ncol=length(w.init))
  thetas <- matrix(0, nrow=niter, ncol=2)
  locs <- matrix(0, nrow=niter, ncol=nrow(clusters))
  
  beta.cur <- beta.init
  w.cur <- w.init
  
  Alist <- lapply(1:nrow(clusters[clusters$mask,]), function(c) {
    inla.spde.make.A(mesh, loc=as.matrix(masterframe[prior.loc$idxlist[[c]],c("east", "north")]))
  })
  print("Begin Burnin")
  for(iter in 1:nburnin) {
    print(iter)
    output <- runone(clusters, masterframe, prior.loc, 
                     beta=beta.cur, w=w.cur, X, Alist, spde, mesh, spcov, mask, idx.masked)
    beta.cur <- output$beta
    w.cur <- output$w
  }
  
  print("End Burnin")
  for(iter in 1:niter) {
    print(iter)
    output <- runone(clusters, masterframe, prior.loc, 
                     beta=beta.cur, w=w.cur, X, Alist, spde, mesh, spcov, mask, idx.masked)
    beta.cur <- betas[iter, ] <- output$beta
    w.cur <- ws[iter, ] <- output$w
    thetas[iter, ] <- output$theta
    locs[iter, ] <- output$loc
  }
  
  list(betas=betas,
       ws=ws,
       thetas=thetas,
       locs=locs)
}


