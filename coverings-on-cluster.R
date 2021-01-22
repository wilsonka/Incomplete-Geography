#### Covering on the cluster

library(maptools)
library(pracma)


### Calcunorthe normalization factors
create.idxlist.disp <- function(masterframe, clusters) {
  distM <- distmat(as.matrix(clusters[,c("east", "north")]), 
                   as.matrix(masterframe[,c("east","north")]))
  
  lapply(1:nrow(clusters), function(x) {
    if(clusters$strata[x] == "R") {
      which(distM[x,] < 10 & masterframe$strata == "R" & 
              masterframe$adminarea == clusters$adminarea[x])
    } else {
      which(distM[x,] < 2 & masterframe$strata == "U" &
              masterframe$adminarea == clusters$adminarea[x])
    }
  })
}




samplept0.5 <- function(pt, adm0.5) {
  r <- runif(1, 0, pt$R)
  theta <- runif(1, 0, 2*pi)
  dx <- r * cos(theta)
  dy <- r * sin(theta)
  new_north  = pt$north + dy
  new_east = pt$east + dx
  newpt <- SpatialPoints(cbind(c(new_east, 0), c(new_north, 0)),
                         proj4string = CRS(proj4string(adm0.5)))
  if(!is.na(over(newpt, adm0.5)$Province[1])  & over(newpt, adm0.5)$Province[1] == pt$adminarea) { # make sure it falls in admin area
    TRUE
  } else {
    FALSE
  }
}

calcDist2ClosestPoly <- function(masterframe, clusters, adm, admtype) {
  if(admtype=="adm0.5") {
    masterframe$adminarea <- masterframe$admin0.5
  } else {
    masterframe$adminarea <- masterframe$admin1
  }
  
  
  idxlist <- create.idxlist.disp(masterframe, clusters)
  
  idx.unique <-sort(unique(unlist(idxlist))) # all possible clusters
  
  tmp <- as(adm, "SpatialLines")
  
  dist2closestpoly <- rep(0, length(idx.unique))
  i <- 1
  for(cutpt in c(seq(1000, length(idx.unique), 1000), length(idx.unique))) {
    print(i)
    pts <- masterframe[idx.unique,c("east", "north")][i:(cutpt),]
    spts <- SpatialPoints(pts, proj4string = CRS(proj4string(adm)))
    dist2closestpoly[i:(cutpt)] <- apply(gDistance(spts, tmp, byid=TRUE),2,min)
    i <- cutpt+1
  }
  
  dist2closestpoly
}

calcNormFactors <- function(masterframe, clusters, admtype, adm, dist2closestpoly,
                            seed=230123L) {
  
  if(admtype=="adm0.5") {
    masterframe$adminarea <- masterframe$admin0.5
  } else {
    masterframe$adminarea <- masterframe$admin1
  }
  
  
  idxlist <- create.idxlist.disp(masterframe, clusters)
  
  idx.unique <- sort(unique(unlist(idxlist))) # all possible clusters
  
  
  masterframe$maxR <- ifelse(masterframe$strata=="U", 2, 10)
  idx.unique.needC <- idx.unique[dist2closestpoly < masterframe$maxR[idx.unique]]
  idx.unique.needC.idx <- which(dist2closestpoly < masterframe$maxR[idx.unique])
  
  set.seed(seed)
  
  covering.needC <- sapply(1:length(idx.unique.needC), function(i) {
    print(i)
    row <- masterframe[idx.unique.needC[i],]
    if(row$strata=="R") {
      # Need to do 2 x
      if(dist2closestpoly[idx.unique.needC.idx[i]] < 5) {
        pt.tmp <- data.frame(east = row$east, north = row$north, R = 5, adminarea = row$adminarea)
        cover5 <- replicate(1000, samplept0.5(pt.tmp, adm))
      } else {
        cover5 <- 1
      }
      
      pt.tmp <- data.frame(east = row$east, north = row$north, R = 10, adminarea = row$adminarea)
      cover10 <- replicate(1000, samplept0.5(pt.tmp, adm))
      c(1/mean(cover5), 1/mean(cover10))
    } else {
      pt.tmp <- data.frame(east = row$east, north = row$north, R = 2, adminarea = row$adminarea)
      cover2 <- replicate(1000, samplept0.5(pt.tmp, adm))
      c(1/mean(cover2), NA)
    }
  })
  
  coverings <- matrix(1, nrow=2, ncol=length(idx.unique))
  coverings[,idx.unique.needC.idx] <- covering.needC
  
  list(coverings=coverings,
       idx.unique=idx.unique)
}



load("Data/forcluster-covering05.rda")
covering05 <- calcNormFactors(masterframe, clusters08, admtype="adm0.5", 
                              adm=adm0.5, dist2closestpoly= dist2closestpoly08,
                              seed=230123L)
save(covering05, file = "Results/covering05.rda")



