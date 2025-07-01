### Create Empty List ###
create_list = function(data = "simulated", rho = c(0.8, 0.5, 0.2), num_train,
                       method =  c("Random", "CD_mean_V2", "CD_mean_MET"),
                       rep = 2000, env = nrow(num_train)){
  
  r = length(rho) ; n = length(num_train) ; met = length(method)
  
  df = sapply(seq(met), function(x){rep(NA, rep)}) %>% data.frame() ; colnames(df) = method
  
  if(data == "simulated"){
    l = lapply(seq(r), function(x){df}) ; names(l) = sprintf("rho = %0.1f", rho)
  } else if(data == "real"){ 
    l = df 
  } else {stop("Error: 'simulated' or 'real'")}
  
  l2 = lapply(seq(n), function(x){l})
  for(i in seq(ncol(num_train))){ names(l2)[i] = paste0("n[",seq(env),  "] = ", num_train[,i], collapse = ", ") }
  
  intersect = l2
  
  index = lapply(seq(env+1), function(x){l2}) ; names(index) = c( paste0("ENV", seq(env)),"Overall")
  datasets = lapply(seq(8), function(x){index}) ; names(datasets) = c("Cor", "CD[mean(v2)]", "CD[mean.MET]", "r^2",
                                                                      "NDCG", "SRC", "RS[ratio]", " ")
  
  return(datasets)
} %>% cmpfun()

###=========###

### MGE Function ###
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

mmer_output = function(K, Nc, output, taxa, env, Trait, g_true, r, iter){
  
  mean_col = function(col) { Reduce(`+`, lapply(Env, function(df) df[[col]])) / length(Env) }
  
  # Population Mean
  b = output[["Beta"]][,3] ; mean_test = b[1]
  for(e in 2:env) {mean_test[e] = b[1]+b[e]}
  
  # Predicted Genetic Effect
  gh = sapply(seq(env+1), function(e){ output[["U"]][[e]][[1]] })[taxa,] %>% data.frame()
  g_pred = sapply(seq(env), function(e){ gh[,1]+gh[,e+1] }) %>% data.frame()
  GEBV = sapply(seq(env), function(e){ g_pred[,e] + mean_test[e]}) %>% data.frame()
  
  
  
  Env = lapply(seq(env), function(e) { 
    df = data.frame(
      TBV = Trait$TBV[which(Trait$Env == e)],
      GEBV = GEBV[, e],
      g = g_true[[r]][Nc * (e - 1) + seq(Nc), iter],
      gh = g_pred[, e],
      Rank_T = rank(-Trait$TBV[which(Trait$Env == e)]),
      Rank_G = rank(-GEBV[, e])
    )
    row.names(df) = row.names(gh)
    return(df)
  })
  
  Env$Overall = data.frame(
    TBV = mean_col("TBV"),
    GEBV = mean_col("GEBV"),
    g = mean_col("g"),
    gh = mean_col("gh"),
    Rank_T = rank(- mean_col("TBV") ),
    Rank_G = rank(- mean_col("GEBV") )
  )
  row.names(Env$Overall) = row.names(gh)
  
  return(Env)
} %>% cmpfun()

mmer_output_real = function(K, Nc, output, taxa, env, Trait, iter){
  
  mean_col = function(col) { Reduce(`+`, lapply(Env, function(df) df[[col]])) / length(Env) }
  
  # Population Mean
  b = output[["Beta"]][,3] ; mean_test = b[1]
  for(e in 2:env) {mean_test[e] = b[1]+b[e]}
  
  # Predicted Genetic Effect
  gh = sapply(seq(env+1), function(e){ output[["U"]][[e]][[1]] })[taxa,] %>% data.frame()
  g_pred = sapply(seq(env), function(e){ gh[,1]+gh[,e+1] }) %>% data.frame()
  GEBV = sapply(seq(env), function(e){ g_pred[,e] + mean_test[e]}) %>% data.frame()
  
  
  
  Env = lapply(seq(env), function(e) { 
    df = data.frame(
      Obs = Trait$Obs[which(Trait$Env == e)],
      GEBV = GEBV[, e],
      Rank_T = rank(-Trait$Obs[which(Trait$Env == e)]),
      Rank_G = rank(-GEBV[, e])
    )
    row.names(df) = row.names(gh)
    return(df)
  })
  
  Env$Overall = data.frame(
    Obs = mean_col("Obs"),
    GEBV = mean_col("GEBV"),
    Rank_T = rank(- mean_col("Obs") ),
    Rank_G = rank(- mean_col("GEBV") )
  )
  row.names(Env$Overall) = row.names(Env[[1]])
  
  return(Env)
} %>% cmpfun()
###=========###

### NDCG ###
ndcg = function(input, top){
  
  # top 是 candidate 數目 * 0.05
  f   = input[order(input$GEBV, decreasing = TRUE), ] %>% .[seq(top), ] # 以 GEBV 排名的前 top 名
  f0  = input[order(input[,1] , decreasing = TRUE), ] %>% .[seq(top), ] # 以 TBV  排名的前 top 名
  d0  = 1/log( seq(top) + 1, base = 2) # 懲罰項公式
  
  NDCG = sum(f[,1]*d0)/sum(f0[,1]*d0)
  
  return(NDCG)
} %>% cmpfun()
###=========###

### Predicting Accuracy ###
simulation = function(K, datasets, g_true, e_true, mean = c(100, 150, 200), iter = seq(2000),
                      sigma_G = 10, sigma_G1 = 3, sigma_G2 = c(2, 20.76, 182), rho = c(0.8, 0.5, 0.2), num_train, DatasetName,
                      method =  c("Random", "CD_mean_V2", "CD_mean_MET"), file_path = "Simulation_Results", 
                      iter_n = seq(ncol(num_train)), iter_r = seq(rho), iter_met = seq(method)){
  
  env = nrow(num_train)
  Nc = nrow(K)
  mu = rep(mean[seq(env)], each = Nc)
  taxa = row.names(K)
  
  # Progress Bar
  suppressPackageStartupMessages({
    library(httr) 
    library(progress)
  })
  
  n_values = apply(num_train, 2, paste, collapse = " & ")

  for (n in iter_n){
    for (r in iter_r){
      cat("\n\n'", "Dataset : ", DatasetName, 
          "    Training Set Size = ", n_values[n], 
          "    rho = ", (rho[r]),
          "'\n", sep = "")
      flush.console()
      
      pb = progress_bar$new(format = paste0(" [:bar] :current/", length(iter), " :percent ; Time: :elapsedfull"), 
                            clear = FALSE, width = 100, total = length(iter))
      
      ### Simulate Training Set 
      for (g in iter){
        used_samples = c()  # 為每個 g 初始化已用的 sample index 記錄
        attempts = 0        # 計數器初始化
        
          repeat {
            if (attempts >= 10) { break }
            
            attempts = attempts + 1
            
            tryCatch({
              sample_idx = sample(setdiff(seq(2000), used_samples), 1)
              used_samples = c(used_samples, sample_idx)
              
              Trait = data.frame(TBV = mu + g_true[[r]][,g] ,
                                 Obs = mu + g_true[[r]][,g] + e_true[[r]][, sample_idx] ,
                                 Env = rep( seq(env), each = Nc ) %>% as.factor() ,
                                 Taxa = rep(taxa, times = env) %>% as.factor())
              
              # 如果 Obs 有 < 0 的數值就重新取樣
              if (any(Trait$Obs < 0)) { next }
              
              for (met in iter_met){
                
                # Select Training Set
                if(met == 1){
                  train_taxa = lapply(seq(env), function(e){ sample(taxa, size = num_train[e, n]) %>% sort() })  
                } else if(met == 2){
                  train_taxa = lapply(seq(env), function(e){ taxa[ Train_CD_mean_v2[[n]][[1]][,e] ] %>% sort() })
                } else if(met == 3){
                  train_taxa = lapply(seq(env), function(e){ taxa[ Train_CD_mean_MET[[n]][[1]][,e] ] %>% sort() })
                }
                
                train_set = lapply(seq(env), function(e){ 
                  Trait[ which(Trait$Taxa %in% train_taxa[[e]] & Trait$Env == e) , ]
                }) %>% do.call(rbind, .)
                train_set$Taxa = as.factor(train_set$Taxa)
                
                MGE = mmer_MGE(K, train_set)
                Env = mmer_output(K, Nc, MGE, taxa, env, Trait, g_true, r, g)
                
                # Indicators
                for (e in seq(env)){ 
                  datasets$`CD[mean(v2)]`[[e]][[n]][[r]][g, met] <- 
                    CD_un_diff(K, train_taxa, index = "CD_mean_V2", var_gxe = c( rep(sigma_G1, (env-1)), sigma_G2[r]), var_g = sigma_G, trial = e) 
                  datasets$`CD[mean.MET]`[[e]][[n]][[r]][g, met] <- 
                    CD_un_diff(K, train_taxa, index = "CD_mean_MET", var_gxe = c( rep(sigma_G1, (env-1)), sigma_G2[r]), var_g = sigma_G, trial = e) 
                }
                datasets$`CD[mean(v2)]`[["Overall"]][[n]][[r]][g, met] <- 
                  CD_un_diff(K, train_taxa, index = "CD_mean_V2", var_gxe = c( rep(sigma_G1, (env-1)), sigma_G2[r]), var_g = sigma_G, trial = 0)
                datasets$`CD[mean.MET]`[["Overall"]][[n]][[r]][g, met] <- 
                  CD_un_diff(K, train_taxa, index = "CD_mean_MET", var_gxe = c( rep(sigma_G1, (env-1)), sigma_G2[r]), var_g = sigma_G, trial = 0)
      
                #intersect[[n]][[r]][g, met] = Reduce(intersect, train_taxa) %>% length() #/ num_train[n]
                
                for (e in seq(env+1)){
                  top = Nc*0.05
                  datasets$`r^2` [[e]][[n]][[r]][g, met] <-  cor(Env[[e]]$g, Env[[e]]$gh)^2
                  datasets$NDCG  [[e]][[n]][[r]][g, met] <-  ndcg(Env[[e]], top)
                  datasets$` `   [[e]][[n]][[r]][g, met] <-  cor(Env[[e]]$TBV, Env[[e]]$GEBV)
                  
                  GEBV_top = Env[[e]][order(Env[[e]]$GEBV, decreasing = TRUE), ] %>% .[seq(top),]
                  datasets$Cor     [[e]][[n]][[r]][g, met] <-  cor(GEBV_top$TBV, GEBV_top$GEBV)
                  datasets$SRC     [[e]][[n]][[r]][g, met] <-  cor(GEBV_top$Rank_T, GEBV_top$Rank_G)
                  # datasets$SRC     [[e]][[n]][[r]][g, met] <-  cor(GEBV_top$TBV, GEBV_top$GEBV, method = "spearman")
                  datasets$`RS[ratio]`[[e]][[n]][[r]][g, met] <-  sum(GEBV_top$Rank_G)/sum(GEBV_top$Rank_T)
                }
              }  
              if ( !all(MGE$U[[1]][[1]] == 0) ) {break}
            },
            error = function(e) { })
          }
    
          pb$tick() 
          save(datasets, file = paste0("數據/", DatasetName, "/", file_path, ".RData"))
        }
        if (!pb$finished) { pb$terminate() }
    }
  }  
  
  return(datasets)
} %>% cmpfun()


real_data = function(K, datasets, Trait, iter = seq(100),
                      num_train, DatasetName,
                      method =  c("Random", "CD_mean_V2", "CD_mean_MET"),
                      iter_n = seq(num_train), iter_met = seq(method), 
                      file_path = "Real_Data_Results"){
  
  env = nrow(num_train)
  Nc = nrow(K)
  taxa = row.names(K)
  
  # Progress Bar
  suppressPackageStartupMessages({
    library(httr) 
    library(progress)
  })
  
  n_values = apply(num_train, 2, paste, collapse = " & ")
  
  for (n in iter_n){
    for (met in iter_met){
      cat("\n\n'", "Dataset : ", DatasetName, 
          "    Training Set Size = ", n_values[n], 
          "    Method = ", method[met], 
          "'\n", sep = "")
      flush.console()
      
      pb = progress_bar$new(format = paste0(" [:bar] :current/", length(iter), " :percent ; Time: :elapsedfull"), 
                            clear = FALSE, width = 100, total = length(iter))
      
      ### Real Data
      if(met == 1){new_iter = iter} else (new_iter = 1)
        
        for (g in new_iter){
          repeat {
            tryCatch({
              
              # Select Training Set
              if(met == 1){
                train_taxa = lapply(seq(env), function(e){ sample(taxa, size = num_train[e, n]) %>% sort() })  
              } else if(met == 2){
                train_taxa = lapply(seq(env), function(e){ taxa[ Train_CD_mean_v2[[n]][[1]][,e] ] %>% sort() })
              } else if(met == 3){
                train_taxa = lapply(seq(env), function(e){ taxa[ Train_CD_mean_MET[[n]][[1]][,e] ] %>% sort() })
              }
              
              train_set = lapply(seq(env), function(e){ 
                Trait[ which(Trait$Taxa %in% train_taxa[[e]] & Trait$Env == e) , ]
              }) %>% do.call(rbind, .)
              train_set$Taxa = as.factor(train_set$Taxa)
              
              MGE = mmer_MGE(K, train_set)
              
              
              Env = mmer_output_real(K, Nc, MGE, taxa, env, Trait, g)
              
              # Indicators
              for (e in seq(env)){ 
                datasets$`CD[mean(v2)]`[[e]][[n]][g, met] <- 
                  CD_un_diff(K, train_taxa, index = "CD_mean_V2", var_gxe = c( rep(10, (env-1)), 10), var_g = 10, trial = e) 
                datasets$`CD[mean.MET]`[[e]][[n]][g, met] <- 
                  CD_un_diff(K, train_taxa, index = "CD_mean_MET", var_gxe = c( rep(10, (env-1)), 10), var_g = 10, trial = e) 
              }
              datasets$`CD[mean(v2)]`[["Overall"]][[n]][g, met] <- 
                CD_un_diff(K, train_taxa, index = "CD_mean_V2", var_gxe = c( rep(10, (env-1)), 10), var_g = 10, trial = 0)
              datasets$`CD[mean.MET]`[["Overall"]][[n]][g, met] <- 
                CD_un_diff(K, train_taxa, index = "CD_mean_MET", var_gxe = c( rep(10, (env-1)), 10), var_g = 10, trial = 0)

              for (e in seq(env+1)){
                top = Nc*0.05
                datasets$`r^2` [[e]][[n]][g, met] <-  cor(Env[[e]]$Obs, Env[[e]]$GEBV)^2
                datasets$NDCG  [[e]][[n]][g, met] <-  ndcg(Env[[e]], top)
                datasets$` `   [[e]][[n]][g, met] <-  cor(Env[[e]]$Obs, Env[[e]]$GEBV)
                
                GEBV_top = Env[[e]][order(Env[[e]]$GEBV, decreasing = TRUE), ] %>% .[seq(top),]
                datasets$Cor     [[e]][[n]][g, met] <-  cor(GEBV_top$Obs, GEBV_top$GEBV)
                datasets$SRC     [[e]][[n]][g, met] <-  cor(GEBV_top$Rank_T, GEBV_top$Rank_G)
                datasets$RS_ratio[[e]][[n]][g, met] <-  sum(GEBV_top$Rank_G)/sum(GEBV_top$Rank_T)
              }
              
            if ( !all(MGE$U[[1]][[1]] == 0) ) {break}
          },
          error = function(e) { })
        }
        
        pb$tick() 
        save(datasets, file = paste0("數據/", DatasetName, "/", file_path, ".RData"))
        }
      if (!pb$finished) { pb$terminate() }
      
    }
  }  
  
  return(datasets)
} %>% cmpfun()
###=========###

### Plot Data ###

plot_data = function(INDEX) {
  list <- list()
  
  for (env_value in names(INDEX)) {
    for (n_value in names(INDEX[[env_value]])) {
      for (rho_value in names(INDEX[[env_value]][[n_value]])) {
        df = INDEX[[env_value]][[n_value]][[rho_value]]
        means = sapply(df, function(x) mean(x, na.rm = TRUE))
        sds = sapply(df, function(x) sd(x, na.rm = TRUE))
        
        data = data.frame(mean = means, sd = sds,
                          env = env_value, n = n_value, 
                          rho = rho_value, Method = names(means))
        list[[paste(env_value, n_value, rho_value)]] = data
      }
    }
  }
  plot_data <- do.call(rbind, list)
  
  
  plot_data$Method <- factor(plot_data$Method, levels = names(means))
  plot_data$env <- factor(plot_data$env)
  plot_data$n <- factor(plot_data$n, levels = names(INDEX[[env_value]]))
  plot_data$rho <- factor(plot_data$rho, levels = rev( names(INDEX[[env_value]][[n_value]]) ))
  
  return(plot_data)
} %>% cmpfun()

create_plot = function(plot_data, title) { 
  
  rho_labels <- setNames(
    sprintf("italic(rho) == %0.1f", levels(data$rho) %>% gsub("rho = ", "", .) %>% as.numeric() ),  ### 有空改一下
    # sprintf("italic(sigma[GX%d]^2) == %.1f", length(levels(plot_data$env))-1,c(10, 20, 30))
    
    levels(plot_data$rho)
  )
  
  n_labels <- setNames(
    gsub("n\\[(\\d+)\\] = (\\d+)", "italic(n)[\\1] == \\2", plot_data$n) %>%
    gsub(",", "~\",\"~", .) ,
    plot_data$n
  )

  # Calculate a single global y-limit based on all data
  overall_max <- max((plot_data$mean+plot_data$sd), na.rm = TRUE)
  overall_min <- min(plot_data$mean, na.rm = TRUE)
  range <- overall_max - overall_min
  y_limits <- c( min(overall_min, overall_min), overall_max )
  
  library(ggplot2) %>% suppressPackageStartupMessages()
  
  ggplot(plot_data, aes(x = env, y = mean, fill = Method)) +
    #ggtitle( parse(text = title) ) +
    labs(subtitle = parse(text = title), size = 50) +
    geom_bar(stat = "identity", alpha = 0.7, position = "dodge", width = 0.8, color = "black") +
    geom_errorbar(aes(ymin = mean, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.8)) +

    # geom_text(aes(label = sprintf("%.2f", mean)),
    #           position = position_dodge(width = 0.8), # Adjust text position to avoid overlap with error bars
    #           vjust = 1.5, # Changed from -0.5 to 1.5 to place text on the bars
    #           size = 3.5) +
    
    facet_grid(rho ~ n, labeller = labeller(rho = as_labeller(rho_labels, default = label_parsed),
                                            n = as_labeller(n_labels, default = label_parsed)), 
               scales = "free_y") + # Added free_y scales, switch = "both"


    coord_cartesian(ylim = y_limits) +  # Add this line to set the y-axis limits for each facet
    
    # scale_color_brewer(palette = "Pastel2") +
    scale_fill_manual(values = c("#9a9a9a", "red", "blue"),
                      labels = c("Random", parse(text = "CD[mean(v2)]"), parse(text = "CD[mean.MET]"))) +  # `V2` 變下標) +
    theme_bw() +

    scale_x_discrete(labels = function(x) {
      ifelse(grepl("^Env\\d+$", x), 
             parse(text = gsub("(Env)(\\d+)", "ENV\\2", x)),
             x)
    }) +
  
    theme(
      plot.title = element_text(size = 60, hjust = 0.5, vjust = -3.0),
      plot.subtitle = element_text(size = 80, hjust = 0.01, vjust = 0.05, face = "bold"),
      plot.margin = margin(60, 20, 20, 20), 
      
      axis.text = element_text(size = 40),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 40),
      axis.title = element_blank(),
      
      strip.text = element_text(size = 40),
      strip.placement = "outside",
      strip.text.x = element_text(size = 25, lineheight = 0.8),
      strip.text.y = element_text(size = 50),
      strip.background = element_rect(color = "black", fill = "lightgrey", linewidth = 1),
      
      legend.title = element_text(size = 60),
      legend.text = element_text(size = 50),
      legend.position = "bottom",
      
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 1.5)  # 加回座標軸的灰框
    )
  
}

add_blank = function(input, option){
  prev_first_val <- NA
  list <- list()
  row = nrow(input) ; col = ncol(input)
  
  if (option == "row"){
    for (i in 1:row) {
      new <- input[i,]
      if (!is.na(prev_first_val) && new[1] != prev_first_val) {
        list[[length(list)+1]] <- matrix(NA, 1, col) %>% as.data.frame()
        colnames(list[[length(list)]]) = colnames(input)
      }
      list[[length(list)+1]] <- new
      prev_first_val <- new[1]
    }
    output <- do.call(rbind, list) %>% as.data.frame()
    colnames(output) = colnames(input)
  }
  
  if (option == "col"){
    for (i in 1:col) {
      new <- input[,i]
      if (!is.na(prev_first_val) && new[1] != prev_first_val) {
        list[[length(list)+1]] <- matrix(NA, row, 1) %>% as.data.frame()
      }
      list[[length(list)+1]] <- new
      prev_first_val <- new[1]
    }
    output <- do.call(cbind, list) %>% as.data.frame()
    colnames(output) = colnames(input)
  }
  
  return(output)
} %>% cmpfun()

###=========###

### Plot Real Data ###
plot_real_data = function(INDEX) {
  list <- list()
  
  for (env_value in names(INDEX)) {
    for (n_value in names(INDEX[[env_value]])) {
        df = INDEX[[env_value]][[n_value]]
        means = sapply(df, function(x) mean(x, na.rm = TRUE))
        
        data = data.frame(mean = means,
                          env = env_value, n = n_value, 
                          Method = names(means))
        list[[paste(env_value, n_value)]] = data
    }
  }
  plot_data <- do.call(rbind, list)
  
  plot_data$Method <- factor(plot_data$Method, levels = names(means))
  plot_data$env <- factor(plot_data$env)
  plot_data$n <- factor(plot_data$n, levels = names(INDEX[[env_value]]))

  return(plot_data)
} %>% cmpfun()

### 分一格 ###
create_real_plot = function(plot_data, title) { 

  shape_value = c(17, 15, 18)
  shape_value = c(shape_value[1: length(levels(plot_data$env))-1], 16)
  
  # Calculate a single global y-limit based on all data
  overall_max <- max(plot_data$mean, na.rm = TRUE)
  overall_min <- min(plot_data$mean, na.rm = TRUE)
  # sd <- mean(plot_data$sd, na.rm = TRUE)
  range <- overall_max - overall_min
  y_limits <- c(overall_min, (overall_max))
  
  library(ggplot2) %>% suppressPackageStartupMessages()
  
  ggplot(plot_data, aes(x = n, y = mean, shape = env, group = interaction(Method, env), color = Method)) +
    ggtitle( parse(text = title) ) +
    labs(subtitle = paste(DatasetName, "dataset"), size = 50) +
    
    geom_line(linewidth = 2) + 
    geom_point(size = 15, na.rm = TRUE) +
    scale_shape_manual(values = shape_value) +

    coord_cartesian(ylim = y_limits) +  # Add this line to set the y-axis limits for each facet

    scale_color_manual(values = c("#9a9a9a", "red", "blue"),
                       labels = c("Random", parse(text = "CD[mean(v2)]"), parse(text = "CD[mean.MET]"))) +
    theme_bw() +
    
    scale_x_discrete(labels = function(x) {
      parse(text = gsub("n\\[(\\d+)\\] = (\\d+)","italic(n)[\\1] == \\2",
        gsub(",", "~\",\"~", x)))}) +
    
    theme(
      plot.title = element_text(size = 60, hjust = 0.5, vjust = -3.0),
      plot.subtitle = element_text(size = 80, hjust = 0.01, vjust = 5, face = "bold"),
      plot.margin = margin(60, 20, 20, 20), 
      
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 30),
      axis.title = element_blank(),
      
      strip.text = element_text(size = 20),
      strip.placement = "outside",
      strip.text.x = element_text(size = 25, lineheight = 0.8),
      strip.text.y = element_text(size = 30),
      strip.background = element_rect(color = "black", fill = "lightgrey", linewidth = 1),
      
      legend.title = element_text(size = 40),
      legend.text = element_text(size = 30),
      legend.position = "bottom",
      
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 1.5)  # 加回座標軸的灰框
    )
}

### 分四格 ###
create_real_plot = function(plot_data, title) { 

  n_labels <- setNames(
    gsub("n\\[(\\d+)\\] = (\\d+)", "italic(n)[\\1] == \\2", plot_data$n) %>%
      gsub(",", "~\",\"~", .) ,
    plot_data$n
  )
  
  # Calculate a single global y-limit based on all data
  overall_max <- max(plot_data$mean, na.rm = TRUE)
  overall_min <- min(plot_data$mean, na.rm = TRUE)
  # sd <- mean(plot_data$sd, na.rm = TRUE)
  range <- overall_max - overall_min
  y_limits <- c(overall_min, overall_max)
  
  library(ggplot2) %>% suppressPackageStartupMessages()
  
  ggplot(plot_data, aes(x = env, y = mean, group = Method, color = Method)) +
    labs(subtitle = parse(text = title), size = 50) +
    geom_line(linewidth = 2) + 
    geom_point(size = 12, na.rm = TRUE) +
    
    facet_grid( ~ n, labeller = labeller(n = as_labeller(n_labels, default = label_parsed)), scales = "free_y") + # Added free_y scales, switch = "both"
    
    coord_cartesian(ylim = y_limits) +  # Add this line to set the y-axis limits for each facet
    
    # scale_color_brewer(palette = "Pastel2") +
    scale_color_manual(values = c("#9a9a9a", "red", "blue"),
                      labels = c("Random", parse(text = "CD[mean(v2)]"), parse(text = "CD[mean.MET]"))) +  # `V2` 變下標) +
    theme_bw() +
    
    scale_x_discrete(labels = function(x) {
      ifelse(grepl("^Env\\d+$", x), 
             parse(text = gsub("(Env)(\\d+)", "ENV\\2", x)),
             x)
    }) +
    
    theme(
      plot.title = element_text(size = 60, hjust = 0.5, vjust = -3.0),
      plot.subtitle = element_text(size = 80, hjust = 0.01, vjust = 0.05, face = "bold"),
      plot.margin = margin(10, 20, 10, 20), 
      
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 45),
      axis.text.y = element_text(size = 30),
      axis.title = element_blank(),
      
      strip.text = element_text(size = 20), #element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(size = 50, lineheight = 0.8),
      strip.text.y = element_text(size = 30),
      strip.background = element_rect(color = "black", fill = "lightgrey", linewidth = 1),
      
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5),
      legend.position = "none",
      
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 1.5)  # 加回座標軸的灰框
    )
}

