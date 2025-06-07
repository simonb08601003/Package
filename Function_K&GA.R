### Filter ###
filterSNP = function(SNP, phe){
  SNP = SNP[which(row.names(SNP) %in% phe$Taxa), ]
  
  missing_percentage <- colSums(is.na(SNP))/nrow(SNP)
  SNP_cleaned <- SNP[, missing_percentage <= 0.10]
  
  column_means <- colMeans(SNP_cleaned, na.rm = TRUE)
  data_filtered <- SNP_cleaned[, column_means <= 0.90 & column_means >= -0.90]
  
  return(data_filtered)
} %>% cmpfun()

### Kinship ###
kinship = function(a, maj = 1, min = -1){  
  library(pbapply)
  
  process_col = function(col) {
    # 如果沒有缺值就略過填補
    if (anyNA(col)) { col[is.na(col)] = maj }
    col[col == maj] = 1
    col[col == mean(c(maj, min))] = 0
    col[col == min] = -1
    return(col)
  }
  a = pbapply(a, 2, process_col)
  
  a <- as.matrix(a) %>% scale()
  K <- (a %*% t(a)) / ncol(a) #tcrossprod
  
  return(K)
} %>% cmpfun()
###==========###

### CD ###
CD_un_diff = function(K, train, index = "CD_mean_V2", var_e = NULL, var_gxe = c( rep(10, length(train)-1), 20 ), var_g = 10, trial = 0){
  
  env = length(train)
  
  # Variance component of MGE model
  omega_g = matrix(var_g, env, env)
  diag(omega_g) = diag(omega_g) + var_gxe
  
  if( is.null(var_e) ){var_e = diag(omega_g)}
  
  # n : number of training set, nc : number of candidates (whole pop.)
  n = c()
  for(e in seq(env)){n[e] = length(train[[e]])}
  nc = nrow(K)
  
  # Gt, Gct, Mt
  Gt = Mt = matrix(0, nrow = sum(n), ncol = sum(n))
  Gct = matrix(0, nrow = nc*env, ncol = sum(n))
  
  for(i in seq(env)){
    idx_i = sum(n[seq(i)]) - n[i] + seq(n[i])
    
    for(j in seq(env)){
      idx_j = sum(n[seq(j)]) - n[j] + seq(n[j])
      idx_all = nc*(i-1) + seq(nc)
      
      Gt[idx_i, idx_j] = omega_g[i,j] * K[train[[i]], train[[j]]] # Gt
      Gct[idx_all, idx_j] = omega_g[i,j] * K[, train[[j]]] # Gct
    }
    
    ni = n[i]
    Mt[idx_i, idx_i] = 1/var_e[i] * (diag(1, ni, ni) - matrix(1/ni, ni, ni)) # Mt
  }
  
  # A and B
  A = Gct %*% solve(Mt %*% Gt + diag(1, sum(n), sum(n))) %*% Mt %*% t(Gct)
  B = omega_g %x% K
  
  if (trial == 0) {
    
    if(index == "CD_mean_V2"){
      
      CD_value = mean(diag(A) / diag(B)) # Diagonal Element of A and B
      
    } else if (index == "CD_mean_MET"){
      
      A_all = B_all = 0
      
      for (i in seq(env)) {
        i_all = nc*(i-1) + seq(nc)
        for (j in seq(env)) {
          j_all = nc*(j-1) + seq(nc)
          
          A_sub = diag(A[i_all, j_all])
          B_sub = diag(B[i_all, j_all])
          
          A_all = A_all + A_sub
          B_all = B_all + B_sub
        }
      }
      CD_value = mean(A_all / B_all)
    }
    
  } else if (trial %in% seq(env)) { ## restriction: number of trials
    
    idx_trial = seq(nc) + nc*(trial-1)
    CD_value = mean( c(diag(A) / diag(B))[idx_trial] )
    
  }
  return(CD_value)
} %>% cmpfun()
###==========###

### OPT-train Untarget Diff ###
opt_un_diff = function(K, n = c(50, 100), option = "Diff", iter = 0, n_iter = 12000, index){  # n : Number of training set
  cat("\n'", "Dataset :", DatasetName, ";",
      "Training Set Size =", paste(n, collapse = " & "), ";", 
      "Index =", index, "'\n")
  
  ## 1. Initialization
  
  N = nrow(K) # Number of All Varieties
  cand = seq(N) # Candidates
  env = length(n) # Number of Trials
  nt = mean(n)
  
  # Length of List (Same or Diff)
  if (option == "Same" & length(unique(n)) == 1) {
    list_len = 1
  } else if (option == "Diff") {
    list_len = seq(env)
  } else {
    stop("Invalid option. Please choose either 'Same' or 'Diff'.")
  }
  
  # Set Solution Size
  if (nt < 50){ ss = 50 } else if (nt >= 50 && nt <= 200){ ss = nt } else (ss = 200)
  
  
  if(iter == 0){
    ##Solution List
    sol = lapply(seq(ss), function(i) {
      sapply(list_len, function(j) {
        sample(c(rep(1, n[j]), rep(0, N - n[j])))
      }) %>% as.data.frame()
    })
    
    max.score = data.frame(matrix(NA, 1, env + 1))
    colnames(max.score) = c("Overall", paste0("Trial", seq_len(env)))
  } else if (iter > 0){      ##### If GA was interrupted, start from where it was paused. #####
    load( paste0("數據/", DatasetName, "/", index, " ", "Number = ", paste(n, collapse = " & "), ".RData") )
  }
  
  
  prev_sol = lapply(seq(ss), function(i) {
    sapply(list_len, function(j) {rep(0, N)}) %>% as.data.frame()
  })
  prev_score = matrix(0, nrow = ss, ncol = env+1)
  
  ## 2. Start Genetic Algorithm
  
  # Progress Bar
  suppressPackageStartupMessages({
    library(httr)
    library(progress)
  })
  
  
  pb <- progress_bar$new(format = paste0(" [:bar] :current/", n_iter, " :percent ; Time: :elapsedfull ; ","Rate: :tick_rate iter/sec", "     :spin"), 
                         clear = FALSE, width = 100, total = n_iter*1.2)
  pb$tick(iter)
  
  # Parallel Computing
  library(parallel)
  
  myCoreNums = detectCores()
  cl <- makeCluster(myCoreNums*0.6) #myCoreNums-1
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages(
    c(library(parallel), library(tidyverse), library(magrittr), 
      library("glmnet"), library("compiler", include.only = "cmpfun")))))
  clusterExport(cl, c("CD_un_diff", "ss", "env", "cand", "K", "N", "n", "list_len", "index"), envir = environment())
  
  ## 3. Calculate Score
  stop = 0
  while(stop == 0){
    iter = iter + 1
    
    # List out the solutions that are changed last iter
    changed_sol = c()
    for(i in seq(ss)){
      if(!(all(sol[[i]] == prev_sol[[i]]))){
        changed_sol = c(changed_sol, i)
      }
    }
    
    
    clusterExport(cl, c("sol", "ss", "changed_sol", "prev_score", "max.score"), envir = environment())
    score = parLapply(cl, seq(ss), function(i) {
      if (i %in% changed_sol) {
        train = lapply(seq(env), function(e) {
          cand[which(sol[[i]][, (e - 1) %% ncol(sol[[i]]) + 1] == 1)]
        })
        res_all = sapply(0:env, function(trial_id) {
          CD_un_diff(K = K, train = train, index = index, trial = trial_id)
        })
        res_all
      } else {
        prev_score[i, ]
      }
    }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    
    colnames(score) = c("Overall", paste0("Trial", 1:env) )
    
    # Check if I didn't calculate the changed solutions
    #if (exists("prev_score")){
    check = c()
    for(i in seq(ss)){ if( all (prev_score[i,]==score[i,])  ){ check = c(check, i) } }
    if( all(setdiff(seq(ss),changed_sol) != check) & length(check)!=0 ){ stop("Didn't calculate the changed solutions.") }
    #}
    
    prev_sol = sol ; prev_score = score    
    max.score[iter, ] = score[which.max(score$Overall), ]
    row.names(max.score)[iter] = iter  
    
    ### Elite and Delete Solution
    elite = which(rank(-score$Overall, ties.method = "max") <= ceiling(ss * 0.1))
    del = which(rank(-score$Overall, ties.method = "min") >= ceiling(ss * 0.6))
    
    ### Print CD
    opt.sol = which(score$Overall == max(score$Overall))
    opt.sol = sol[[ opt.sol[1] ]]
    opt = lapply(seq(env), function(e) {
      cand[which(opt.sol[, (e - 1) %% ncol(sol[[i]]) + 1] == 1)]
    })
    
    clusterExport(cl, c("opt"), envir = environment())
    max.row = which.max(score$Overall)
    max.score[iter, ] = score[max.row, ]
    
    ### Stop criteria
    th = 0.0001
    if(iter > n_iter){ 
      if (all(sapply(env+1, function(e) { 
        score[max.row, e] - max.score[iter - 500, e] < th}
      ) )) { stop = stop + 1 } 
    } # Converge
    
    if(iter > 2){
      if (all(sapply(1, function(e) { 
        score[max.row, e] < max.score[iter-1, e]}
      ) ))  {stop("Max score drops!")} 
    }
    
    
    
    ## 4. Crossover
    for(i in del){
      repeat {
        chr = sample(seq(ss)[-del], 2)
        pos = sample(seq(N - 1), 1)
        for(e in seq(list_len)){sol[[i]][,e] = c(sol[[ chr[1] ]][1:pos, e], sol[[ chr[2] ]][ (pos+1):N, e])}
        
        # In case n.sol is wrong during crossover
        n.sol = sapply(seq(env), function(e){ length(which(sol[[i]][, (e - 1) %% ncol(sol[[i]]) + 1] == 1)) })
        
        if (all(abs(n.sol-n) == 0)) { break }
      }
    }
    
    ## 5. Mutation    
    sol.new = sol
    r = 0.04*n
    for(i in seq(ss)){
      for(e in seq(list_len)){
        pos = c( sample(which(sol[[i]][,e] == 0), ceiling(r[e])), sample(which(sol[[i]][,e] == 1), ceiling(r[e])) )
        sol.new[[i]][pos, e] = abs(sol.new[[i]][pos, e] - 1)
      }
    }
    
    clusterExport(cl, c("ss", "sol.new", "sol", "elite", "del", "score", "n"), envir = environment())   
    sol = parLapply(cl, seq(ss), function(i){
      
      if (!(i %in% elite)) {
        sol.new[[i]]  # Return the mutated solution for non-elite cases
        
      } else {
        old = score[i, ]
        train.new = lapply(seq(env), function(e) {
          cand[which(sol.new[[i]][, (e - 1) %% ncol(sol[[i]]) + 1] == 1)]
        })
        
        new = sapply(0:env, function(e) {
          CD_un_diff(K = K, train = train.new, index = index, trial = e)
        })
        
        if ( all(new > old) ) {sol.new[[i]]} else {sol[[i]]}
        
      }
    })
    
    pb$tick()  # Update Progress Bar
    
    save(iter, sol, max.score,
         file = paste0("數據/", DatasetName, "/", index, " ", "Number = ", paste(n, collapse = " & "), ".RData"))
    
    if(stop != 0) {cat("\nGenetic Algorithm has ended!\n\n")}
  }### End GA
  
  stopCluster(cl)
  if (!pb$finished) { pb$terminate() }
  
  ## Output Optimized Training set
  opt.sol = which(score$Overall == max(score$Overall))
  opt.sol = sol[[ opt.sol[1] ]]
  opt.set = sapply(seq(env), function(e) {
    c(cand[which(opt.sol[, (e - 1) %% ncol(sol[[i]]) + 1] == 1)], rep(NA, max(n) - n[e]))
  }) %>% as.data.frame()
  
  return(list(Opt.set = opt.set, CD = max.score))
  
} %>% cmpfun()
###=========###

### Sample Random Effects ###
RandomEffects = function(K, env = 2, OmegaG_fixed = 10, OmegaG_first = 10, OmegaG_last = c(0, 10, 20), 
                         n = 2000, DatasetName, core = 0.6, file_path = "/RandomEffects.RData"){
  
  g_true = e_true = list()
  
  suppressPackageStartupMessages({
    library(MASS, include.only = 'mvrnorm')
    library(parallel)
    library(pbapply)
  })
  
  # Parallel Computing
  myCoreNums = detectCores()
  cl <- makeCluster(myCoreNums*core)
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages(
    c(library(parallel), library(tidyverse), library(magrittr), 
      library("glmnet"), library("compiler", include.only = "cmpfun"),
      library(MASS), library(pbapply)))))
  
  for (r in seq(OmegaG_last) ){
    cat("### ", DatasetName, " ; OmegaG_last = ", (OmegaG_last[r]+OmegaG_fixed[1]), " ###\n", sep = "")
    
    omega_g = matrix(OmegaG_fixed, env, env) + diag(c( rep(OmegaG_first, env-1) , OmegaG_last[r]))
    omega_e = diag(OmegaG_fixed, env) + diag(c( rep(OmegaG_first, env-1) , OmegaG_last[r]))
    
    G = omega_g %x% K
    E = omega_e %x% diag(nrow(K))
    
    clusterExport(cl, c("G", "E"), envir = environment())
    
    cat("    Simulated Genetic Effects\n")
    g_true[[r]] = pbsapply(seq(n), cl = cl, function(i) {
      mvrnorm(n = 1, mu = rep(0, nrow(G)), Sigma = G)
    }) %>% as.data.frame()
    save(g_true, e_true, file = paste0("數據/", DatasetName, file_path))
    
    cat("    Simulated Environmental Effects\n")
    e_true[[r]] = pbsapply(seq(n), cl = cl, function(i) {
      mvrnorm(n = 1, mu = rep(0, nrow(E)), Sigma = E)
    }) %>% as.data.frame()
    save(g_true, e_true, file = paste0("數據/", DatasetName, file_path))
  }
  stopCluster(cl)
  return(list(g_true, e_true))
  
} %>% cmpfun()
###=========###

change = function(set, order = c(3,2,1), n){
  new = set
  for(i in n){
    for(o in seq(order)){
      new[[i]]$Opt.set[, o] = set[[i]]$Opt.set[, order[o]]
      new[[i]]$CD[, o] = set[[i]]$CD[, order[o]]
    }
  }
  return(new)
}

### Original Version ###
# kinship_alt = function(a){
#   
#   # Impute missing data as Major allele
#   for (i in 1:length(a)){
#     a[which(is.na(a[,i])),i]  = 1
#     a[which(a[,i] == -1), i]  = 0
#   }
#   
#   # Calculate MAF
#   p = apply(a, 2, function(x) mean(x,na.rm = F)) %>% data.frame() %>% t()
#   
#   # Number of individuals
#   ind = length(c(a[,1]))
#   
#   # n - 2p
#   for (i in 1:length(a)){
#     a[which(a[,i] == 1),i]    = 2 - 2 * p[i]
#     a[which(a[,i] == 0),i]    = 0 - 2 * p[i]
#   }
#   a = as.matrix(a)
#   K = a %*% t(a)  / (tr(a %*% t(a))/ind)
#   K = as.matrix(K)
#   
#   return(K)
# } %>% cmpfun()


# CD_un_diff_alt = function(K, train, index = "CD_mean_V2", h2 = 0.5, var_gxe = c( rep(1, length(train)-1), 0.5 ), cov_g = 0.5, trial = 0){
#   
#   env = length(train)
#   
#   # Variance component of MGE model
#   omega_g = matrix(cov_g, env, env) + diag(var_gxe)
#   var_g = var_gxe + cov_g
#   var_e = (1 - h2)/h2 * var_g
#   
#   
#   # nt : number of training set, nc : number of candidates (whole pop.)
#   n = c()
#   for(e in seq(env)){n[e] = length(train[[e]])}
#   nc = nrow(K)
#   
#   # Gt
#   Gt = matrix(0, nrow = sum(n), ncol = sum(n))
#   
#   for(i in seq(env)){
#     for(j in seq(env)){Gt[sum(n[seq(i)]) - n[i] + seq(n[i]), sum(n[seq(j)]) - n[j] + seq(n[j])] = cov_g* K[train[[i]], train[[j]]]}
#   }
#   
#   
#   for(e in seq(env)){Gt[sum(n[seq(e)]) - n[e] + seq(n[e]), sum(n[seq(e)]) - n[e] + seq(n[e])] = var_g[e]* K[train[[e]], train[[e]]]}
#   
#   # Gct
#   Gct = matrix(0, nrow = nc*env, ncol = sum(n))
#   
#   for(e in seq(env)){
#     for(j in seq(env)){
#       Gct[nc*(e-1) + seq(nc), sum(n[seq(j)]) - n[j] + seq(n[j])] = cov_g * K[, train[[j]]]
#     }
#   }
#   for(e in seq(env)){Gct[nc*(e-1) + seq(nc), sum(n[seq(e)]) - n[e] + seq(n[e])] = var_g[e] * K[, train[[e]]]}
#   
#   # IJ: I - J_bar
#   IJ = function(x){diag(1, x, x) - matrix(1/x, x, x)}
#   Mt = matrix(0, nrow = sum(n), ncol = sum(n))
#   for(e in seq(env)){Mt[sum(n[seq(e)]) - n[e] + seq(n[e]), sum(n[seq(e)]) - n[e] + seq(n[e])] = 1/var_e[e]* IJ(n[e])}
#   
#   # Diagonal Element of A and B
#   A = Gct %*% solve(Mt %*% Gt + diag(1, sum(n), sum(n))) %*% Mt %*% t(Gct)
#   B = omega_g %x% K
#   
#   if (trial == 0) {
#     
#     if(index == "CD_mean_V2"){
#       CD_value = mean(diag(A) / diag(B))
#     } else if (index == "CD_average"){
#       CD_value = sum(A) / sum(B)
#     }
#     
#   } else if (trial %in% seq(length(var_g))) { ## restriction: number of trials
#     
#     if(index == "CD_mean_V2"){
#       CD_value = mean( c(diag(A) / diag(B))[seq(nc) + nc*(trial-1)] )
#     } else if (index == "CD_average"){
#       CD_value = mean(A[seq(nc) + nc*(trial-1) ,seq(nc) + nc*(trial-1) ]) / mean(B[seq(nc) + nc*(trial-1) , seq(nc) + nc*(trial-1)])
#     }
#   }
#   return(CD_value)
# } %>% cmpfun()