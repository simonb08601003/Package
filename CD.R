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
