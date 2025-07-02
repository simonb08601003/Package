GA = function(K, n = c(50, 100), option = "Diff", iter = 0, n_iter = 12000, index){  # n : Number of training set
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
