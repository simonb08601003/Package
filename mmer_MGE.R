mmer_MGE = function(K, train_set){
  #train_set: TBV, Obs, Env, Taxa
  library(sommer) %>% suppressPackageStartupMessages()
  
  output = mmer(
    Obs ~ Env,  # Fixed Effect
    random = ~ vsr(Taxa, Gu = K) + vsr(dsr(Env), Taxa, Gu = K),  # Random Effect
    rcov = ~ vsr(dsr(Env), units),  # Residual Matrix
    data = train_set,  # Training Set
    verbose = FALSE,
    nIters = 10000,
    method = "AI",
    dateWarning = F
  )
  return(output)
} %>% cmpfun()
