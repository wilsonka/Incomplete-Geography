#### Setup prior for cluster


create.priorlist.disp <- function(masterframe, clusters, admtype="adm0.5", covering = NULL) {
  ## DISPLACEMENT
  ##    by long, lat contains N_{EA} * 1/(2 pi R d) where d is the distance b/w
  ##        long, lat and EA, R is 2 for urban and 5 for rural.
  ##    NOTE: For rural, 1% are further displaced up to 10km
  
  if(admtype=="adm0.5") {
    masterframe$adminarea <- masterframe$admin0.5
  } else {
    masterframe$adminarea <- masterframe$admin1
  }
  
  idxlistlist <- create.idxlist.disp(masterframe, clusters, returndist = T)
  idxlist <- idxlistlist$idxlist
  distlist <- idxlistlist$distlist
  
  coverings <- covering$coverings
  idx.unique <- covering$idx.unique
  
  
  prior <- lapply(1:length(idxlist), function(x) {
    id <- idxlist[[x]]
    if(length(id) == 1) {
      1
    } else {
      
      Nisk <- 1 # Change if the EAs have different sizes that influence sampling
      
      dists <- distlist[[x]]
      
      if(!is.null(covering)) {
        id_to_idxunique <- sapply(1:length(id) , function(i) which(idx.unique == id[i]))
        cs <- coverings[,id_to_idxunique]
      } else {
        cs <- matrix(1, nrow=2, ncol=length(id))
      }
      if(clusters$strata[x] == "R") {
        unnormalized <- Nisk * (0.99 * (dists < 5) * cs[1,] * 1/(2*pi * 5 * dists) +
                                  0.01 * cs[2, ] * 1/(2 * pi * 10 * dists))
        
      } else {
        unnormalized <- Nisk * cs[1, ] * 1/(2*pi * 2 * dists)
      }
      unnormalized/sum(unnormalized)
    }
    
  })
  return(list(idxlist=idxlist,
              prior=prior))
}


create.idxlist.mask <- function(masterframe, clusters) {
  idxlist <- lapply(1:nrow(clusters), function(x) {
    which(masterframe$adminarea == clusters$adminarea[x] &
            masterframe$strata == clusters$strata[x])
  })
  return(idxlist)
}




create.priorlist.mask <- function(masterframe, clusters, admtype="adm0.5") {
  ## MASKING
  ##    by long, lat contains N_{EA}
  
  if(admtype=="adm0.5") {
    masterframe$adminarea <- masterframe$admin0.5
  } else {
    masterframe$adminarea <- masterframe$admin1
  }
  
  idxlist <- create.idxlist.mask(masterframe, clusters)
  
  
  
  prior <- lapply(1:length(idxlist), function(x) {
    id <- idxlist[[x]]
    if(length(id) == 1) {
      1
    } else {
      
      Nisk <- 1 # Change if the EAs have different sizes that influence sampling
      unnormalized <- rep(Nisk, length(id))
      
      unnormalized/sum(unnormalized)
    }
    
  })
  return(list(idxlist=idxlist,
              prior=prior))
}
