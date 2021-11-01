# This is a filter that needs the custom geodiff package in this branch to run. 
# Basically it tries to fit the data in each fov to the NBthDE model
# for each fov that has columns failing more than a threshold (10% here) 
# that fov is excluded from the subsequent analyses

filter_fovs = function(annot,
                       raw0,
                       features_high,
                       features_all,
                       sizefact,
                       sizefact0,
                       gamma0){
  fov_start = min(unique(annot$fov))
  fov_end = max(unique(annot$fov))
  suitable_fovs <- vector(length = fov_end - fov_start +1)
  failure_threshold = 10
  for (x in seq(fov_start,fov_end)){
    print(x)
    high_ROIs <- names(which(sizefact[annot$fov%in%c(x)]>0))
    
    
    
    
    features_all <- rownames(raw0)
    
    
    
    
    NBthDEmod2 <- fitNBthDE(form = ~fov,
                            annot=annot[high_ROIs, ],
                            object=raw0[, high_ROIs],
                            probenum = rep(1, NROW(raw0[, high_ROIs])),
                            features_high = features_high,
                            features_all = features_all,
                            sizefact_start=sizefact[high_ROIs],
                            sizefact_BG=sizefact0[high_ROIs],
                            threshold_mean = gamma0,
                            preci2=10000,
                            prior_type="contrast",
                            covrob=FALSE,
                            preci1con=1/25,
                            sizefactrec=FALSE,
                            sizescalebythreshold=TRUE)
    percentage_failed_columns = 100*sum(is.na(NBthDEmod2$conv))/length(NBthDEmod2$conv)
    print("percentage of failed columns:")
    print(percentage_failed_columns)
    if (percentage_failed_columns<failure_threshold){
      suitable_fovs[x-fov_start+1] = x
    }
  }
  suitable_fovs = suitable_fovs[suitable_fovs>0]
  return(suitable_fovs)
}
