# this module takes accepts the GeoDiff
source("~/Downloads/Seurat_GeoDiff_object.R")
source("~/r-notebooks/nanostring_cpp_opti/just_rcpp/R/fov_filter.R")

GeoDiff_DE<- function(Seurat_GeoDiff_Obj){
  
  neg<- Seurat_GeoDiff_Obj[["neg"]]
  # extracting the negative probes from the object
  
  raw<- Seurat_GeoDiff_Obj[["raw"]]
  # extracting the raw counts from the object
  
  annot <- Seurat_GeoDiff_Obj[["annot"]]
  # extracting the annotation from the object
  
   
  mu <- median(rowSums(neg))
  # For each negative probe, calculate total count for all cells.Calculate the median mu.
  
  
  pos_count <- rowSums(raw)
  # For each positive probe, calculate total count for all cells.
  
  
  indx_low_pos <- which(pos_count < mu)
  #Select positive probes with total count less than mu, call them low positive probes
  
  negmod <- GeoDiff:::fitPoisBG(rbind(neg, raw[indx_low_pos, ]), size_scale = "sum") # use sum in SMI data. use first will lead to distortion in data
  negdiag2 <- GeoDiff:::diagPoisBG(negmod, generate_ppplot = FALSE)
  # Combine low positive probes and negative probes, fit the Poisson Background model, implement the common diagnostics procedure
  
  
  negmod2 <- GeoDiff:::fitPoisBG(neg, size_scale = "sum") # use sum in SMI data. use first will lead to distortion in data
  negmod2$sizefact <- negmod$sizefact
  
  
  sc <- GeoDiff:::BGScoreTest(raw, negmod2, adj = 1, removeoutlier = FALSE, useprior = TRUE)
  # perform score tests
  
  
  high_features <- ((sc$scores > quantile(sc$scores, probs = 0.4)) & (sc$scores < quantile(sc$scores, probs = 0.95))) %>%
    which() %>%
    names()
  # Selecting high features
  
  sizefact0 <- negmod$sizefact
  gamma0 <- mean(negmod2$featfact)
  gamma_features <- rowSums(raw[high_features,]) - gamma0
  sizefact <- (colSums(raw[high_features,])-length(high_features)*gamma0*sizefact0)/sum(gamma_features)
  # calculate the sizefact
  # estimate a_j from the poisson threshold model
  
  
  sizefact[sizefact<0] <- 0
  # replace negative by 0
  
  
  gamma_features <- gamma_features*sum(sizefact)
  sizefact <- sizefact/sum(sizefact)
  # rescale
  
  ################## DE #################
  
  rownames(annot) <- colnames(raw)
  # changing the row names of the annotation file to correspond to the colnames of the raw count
  
  #high_ROIs <- names(which(sizefact[annot$fov%in%c(1:length(unique(unlist(annot$fov))))]>0))
  
  suitable_fovs = filter_fovs(annot,
                             raw,
                             high_features,
                             features_all,
                             sizefact,
                             sizefact0,
                             gamma0)
  
  high_ROIs <- names(which(sizefact[annot$fov%in%c(suitable_fovs)]>0))
  # extracting cells based on the sizefact
  
  
  
  features_all <- rownames(raw)
  # extracting the row names for the raw count
  
  NBthDEmod2 <- GeoDiff:::fitNBthDE(form = ~fov,
                                    annot=annot[high_ROIs, ],
                                    object=raw[, high_ROIs],
                                    features_high = high_features,
                                    features_all = features_all,
                                    sizefact_start=sizefact[high_ROIs],
                                    sizefact_BG=sizefact0[high_ROIs],
                                    probenum = rep(1, NROW(raw[, high_ROIs])),
                                    threshold_mean = gamma0,
                                    preci2=10000,
                                    prior_type="contrast",
                                    covrob=FALSE,
                                    preci1con=1/25,
                                    sizefactrec=FALSE,
                                    sizescalebythreshold=TRUE)
  
  # fitting the data using the fitNBthDE
  
  return(NBthDEmod2)
  # returning the output of the DE model 
  
}
