rm(list=ls())
set.seed(123)  # for reproducibility
library(Matrix)
library(MASS)
library(glmnet)        # Elastic Net
library(ncvreg)        # SCAD
library(monomvn)       # BLASSO
library(MBSP)          # MBSP
library(BayesSUR)      # BayesSUR
library(LaplacesDemon) # sampling omega

tau <- 1e-5
# Helper Functions 
construct_x <- function(p, type,strong_rows) {
  if (type == "indpt") {
    X = matrix(rnorm(n*p),n,p)
    X = scale(X)
    return(X)
  } else if (type == "AR1") {
    cov_x <- generate_cov_x_ar1(p, rho = 0.7)
    X <- mvrnorm(n = n, mu = rep(0, p), Sigma = cov_x)
    X <- scale(X, center=TRUE, scale=TRUE)
    return(X)
  } else if (type == "block") {
    cov_x <- generate_cov_x_block(p, rho = 0.7)
    X <- mvrnorm(n = n, mu = rep(0, p), Sigma = cov_x)
    X <- scale(X, center=TRUE, scale=TRUE)
    return(X)
  } else if (type == "diag") {
    cov_x <- generate_cov_x_diag(p, p0 = p0, strong_rows)
    X <- mvrnorm(n = n, mu = rep(0, p), Sigma = cov_x)
    X <- scale(X, center=TRUE, scale=TRUE)
    return(X)}
}

generate_cov_x_ar1 <- function(p, rho) {
  times <- 1:p
  H <- abs(outer(times, times, "-"))
  cov_x <- rho^H
  return(cov_x)
}

generate_cov_x_diag <- function(p,p0,indx.beta=strong_rows){
    Z1 <- diag(p0)
    Z <- diag(p - p0)
    Z2 <- 0.25 * (Z + p0 * rep(1, p - p0)%*%t(rep(1, p - p0)))
    # library(Matrix)
    final_mat <- matrix(0, p, p)
    
    # Relevant block: identity on indx.beta
    final_mat[indx.beta, indx.beta] <- Z1
    
    all_idx <- 1:p
    irr_idx <- setdiff(all_idx, indx.beta)
    final_mat[irr_idx, irr_idx] <- Z2
    return(final_mat)
}


generate_cov_x_block <- function(p, rho) {
  block_sizes <- c(30, 40, 20, 110, 200, 200, 200, p - 800)
  final_mat <- matrix(rho, nrow = block_sizes[1], ncol = block_sizes[1]) + diag(1 - rho, block_sizes[1])
  for (i in block_sizes[-1]){
    a <-  matrix(rho, nrow = i, ncol = i) + diag(1 - rho, i)
    final_mat <- (as.matrix(bdiag(final_mat, a)))
  }
  
  # Stabilize if needed
  eig_vals <- eigen(final_mat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig_vals <= 0)) diag(final_mat) <- diag(final_mat) + 1e-6
  
  return(final_mat)
}

construct_B <- function(p, q, p0) {
  B <- matrix(0, nrow = p, ncol = q)
  strong_rows <- sort(sample(1:p, p0))
  top <- strong_rows[1:floor(p0/3)]
  mid <- strong_rows[ceiling(p0/3):floor(2*p0/3)]
  low <- strong_rows[ceiling(2*p0/3):p0]
  
  B[top, ] <- matrix(sample(c(-2, -1, 1, 2), length(top)*q, prob=c(0.1,0.4,0.4,0.1),replace = TRUE), nrow = length(top))
  
  for (r in mid) {
    pos <- sample(1:q, q/2)
    B[r, pos] <- sample(c(-1,-2, 2,1),prob=c(0.2,0.3,0.3,0.2), length(pos), replace = TRUE)
  }
  
  for (r in low) {
    pos <- sample(1:q, ceiling(q/4))
    B[r, pos] <- sample(c(-2, 2), length(pos), replace = TRUE)
  }
  norms_q <- NULL
  for ( i in strong_rows){
    norms_q <- rbind(norms_q, norm(B[i,1:q],type = '2'))
  }
  B1 <- B/floor(sqrt(max(norms_q)))
  
  return(list(B=B1,
              indx.beta=strong_rows))
}

generate_cov_y <- function(q, type) {
  if (type == "dense") {
    rho <- 0.7
  } else {
    rho <- 0.15
  }
  return(rho * matrix(1, q, q) + (1 - rho) * diag(q))
}

# to check B_0
non_zero_counts <- apply(B, 1, function(row) sum(row != 0))
non_zero_counts[indx.beta]

###############################################
# setups 

cov_setups <- list(
  list(p=1000,p0=15,q=10,n=200,cov_x = "indpt", cov_y = "dense"),
  list(p=1500,p0=20,q=10,n=200,cov_x = "AR1",   cov_y = "light"),
  list(p=1000,p0=15,q=10,n=200,cov_x = "block", cov_y = "dense"),
  list(p=1500,p0=12,q=8,n=100,cov_x = "diag", cov_y = "light")
)

# all combinations
#library(MASS)

###############3#
#DO NOT RUN AGAIN!!!!!!!!!!!!!!
for (setup_id in 1:4) {
  s <- cov_setups[[setup_id]]
  n <- s$n
  p <- s$p
  p0 <- s$p0
  q <- s$q
  
  
  # === Simulate X and B ===
  B_gen <- construct_B(p, q, p0)
  B <- B_gen$B
  strong_rows <- B_gen$indx.beta
  X <- construct_x(p, s$cov_x, strong_rows)
  
  
  # === Generate E and Y ===
  cov_y <- generate_cov_y(q, s$cov_y)
  
  ### random omega
  #V <- diag(q)
  #library(MCMCpack)
  #Omega <- riwish(nu, V)
  #Z <- matrix(rnorm(n * q), n, q)
  #E <- Z %*% chol(Omega) 
  
  
  E <- mvrnorm(n, rep(0, q), cov_y)
  Y <- scale(X %*% B + E)
  # subgaussian errors
  #y = cros + rmt(n,S=diag(q)*sqrt(1.5),df=10)
  #y = cros + runif(n,-1,1)
  #y = cros + rbinom(n,1,0.5)*diag(1.5,nrow=q, ncol= q)
  
  
  # === Save files with index ===
  prefix <- paste0( p0, "_", q, "_", setup_id)
  #write.csv(X, paste0("X_", prefix, ".csv"), row.names = FALSE)
  #write.csv(B, paste0("B_", prefix, ".csv"), row.names = FALSE)
  #write.csv(Y, paste0("Y_", prefix, ".csv"), row.names = FALSE)
  #write.csv(strong_rows, paste0("S_", prefix, ".csv"), row.names = FALSE)
}




##################################################################
#Our functions
### log posterior function
log_post <- function(selected_vars){
  pgamma <- length(selected_vars)
  V <- diag(q)
  c <- q + 1
  
  a <- 1/(1 + 1/g)
  b <- 1 + 1/tau
  x <- X[, selected_vars, drop = FALSE]
  P <- x %*% solve(t(x) %*% x) %*% t(x)
  zai <- (1 - 1/b) * diag(n) - (b + 1/b - 2) * a * P / (b - a)
  
  log_num <- (pgamma * q / 2) * log(1 - a) + 
    (c / 2) * log(det(V)) +
    (n * q / 2) * log(1 - 1/b)
  
  log_den <- log(choose(p, pgamma)) +
    (pgamma * q / 2) * log(1 - a / b) + 
    ((n + c) / 2) * log(det(V + t(Y) %*% zai %*% Y))
  
  return(log_num - log_den)
}

### log posterior for comparison
log_post_comparison <- function(selected_vars, s){
  pgamma <- length(selected_vars)
  V <- diag(q)
  c <- q + 1
  a <- 1/(1 + 1/g)
  b <- 1 + 1/tau
  x <- X[, selected_vars, drop = FALSE]
  xc <- X[, c(selected_vars, s), drop = FALSE]
  P <- x %*% solve(t(x) %*% x) %*% t(x)
  Pc <- xc %*% solve(t(xc) %*% xc) %*% t(xc)
  w <- (b + 1/b - 2) * a / (b - a)
  zai <- (1 - 1/b) * diag(n) - w * P
  term <- diag(q) - w * solve(V + t(Y) %*% zai %*% Y) %*% (t(Y) %*% (Pc - P) %*% Y)
  
  log_num <- log(pgamma + 1) + (q / 2) * log(1 - a)
  log_den <- log(p - pgamma) + (q / 2) * log(1 - a / b) + ((n + c) / 2) * log(det(term))
  
  return(log_num - log_den)
}

### post processing functions
PRFJ <- function(model, true){
  TP <-	length(intersect(model, true))#correctly selected	
  FP <- length(setdiff(model, true))#wrongly selected	
  FN <-	length(setdiff(true, model))#missed true vars
  
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- 2*(precision*recall)/(precision+recall)
  Jaccard <- length(intersect(model,true))/
    length(union(model,true))
  norm_BbyB0 <- norm(B[model,],"2")/norm(B[true,],"2")
  
  return(list(precision = precision,
              recall = recall,
              F1 = F1,
              Jaccard_index = Jaccard,
              norm_BbyB0 = norm_BbyB0
              ))
}

psi <-function(x){
  a <- (exp(x) - 1)/(exp(x) + 1)
  return(a)
  }

logi <- function(x) {
  1 - exp(-1000 * x)
}#working fn for alpha

log_scaled <- function(x) {
  vi <- NULL
  for (i in x){
  x <- max(i, 1e-100)  # avoid log(0)
  val <- log(x)
  vi <- cbind(vi,(val - log(1e-100)) / ( - log(1e-100)))  # scale to [0,1]
  }
  return(vi)
  }#not being used


G_gamma <- function(selected_vars){##MN
  a <- 1/(1+ 1/g)
  b <- 1+ 1/tau
  x <- X[, selected_vars]
  P <- x%*%solve(t(x)%*%x)%*%t(x)
  #posterior mean for G_gamma_bar
  G <- kronecker(diag(q),solve(b*diag(n)-a*P)%*%(diag(n)-a*P))%*%as.vector(Y)
  me <- kronecker(diag(q),solve(t(x)%*%x)%*%t(x))
  mean_mat <- me%*%(as.vector(Y)-G)
  mean_mat <- matrix(mean_mat,ncol = q)
  return(mean_mat)
}

library(matrixsampling)# or LaplacesDemon to sample MxVT
B_gamma <- function(selected_vars){
  ### from MxVT(M,Sigma,Omega,neu)
  a <- 1/(1+ 1/g)
  b <- 1+ 1/tau
  c <- q + 1
  pgamma <- length(selected_vars)
  v_n <- 1-b^{-1}
  u_n <- a^{1} - b^{-1}
  x <- X[, selected_vars]
  P <- x%*%solve(t(x)%*%x)%*%t(x)
  # posterior mean for B_gamma
  U <- solve(t(x)%*%x)/u_n
  M <- v_n*U%*%t(x)%*%Y
  V <- diag(q)
  V <- {v_n*t(Y)%*%(diag(n)-(v_n/u_n)*P)%*%Y +V}/(n+c-q+1)
  V_fixed <- as.matrix(nearPD(V)$mat)
  library(matrixsampling)
  B_g <- rmatrixt(n=10, nu=(n+c-q+1), M=M, U=U, V=V_fixed, checkSymmetry = FALSE, keep = TRUE)
  B <-  apply(B_g, c(1, 2), mean)
  return(B)
}
### shouldn't we change the name?
s_53_log <- function(X, Y, temp = seq(1, 0.1, length.out = 20) , 
                     L = 20,J=20, burnin = 0.5,g = g, M =n/log(n)){
  
  #scaling
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  tau <- 1e-5
  
  
  # Initialize empty set of selected variables
  selected_vars <- integer(0)
  model_a_store <- NULL
  model_r_store <- NULL
  accepted_models <- data.frame()
  
  # neighborhood function
  get_neighborhood <- function(selected_vars, X, Y) {
    model <- lm(Y ~ X[, selected_vars, drop=FALSE])  # Fit the model
    r <- residuals(model)
    h <- t(r) %*% X
    h_norms <- apply(h, 2, function(col) norm(col, type = '2'))
    scrs_vars <- order(h_norms, decreasing = TRUE)[1:M]
    setdiff(scrs_vars, selected_vars)  # Return only unselected variables
  }
  
  logsumexp <- function(log_vals) {
    m = max(log_vals)
    return(m + log(sum(exp(log_vals - m))))
  }
  
  
  tic <- proc.time()[3]
  
  # choosing the initial model
  model <- lm(Y ~ 1) # k(1,1)
  r <- residuals(model)
  h <- t(r)%*%X
  h_norms <- apply(h, 2,function(col) norm(col,type = '2'))
  scrs_vars <- order(h_norms, decreasing = TRUE)[1:M]
  remaining_vars <- setdiff(scrs_vars, selected_vars)
  
  # Only consider additive neighborhood for the initial model
  for (s in remaining_vars) {
    p_addition <- log_post(s)  # Compute the marginal posterior probability
    model_a_store <- rbind(model_a_store, data.frame(p_addition,s))
  }
  selected_vars <- model_a_store$s[which.max(model_a_store$p_addition)]#map
  accepted_models <- data.frame(post = max(model_a_store$p_addition),
                                model = selected_vars,
                                add = NA,
                                rem = NA)
  iter <- 0
  iter.before.burn <- 0
  
  foo <- ceiling((L*J)*burnin)
  #starting the loop
  for (l in 1:L) {
    for (i in 1:J){
      iter <- (l - 1) * J + i
      print(paste(l,i,iter,iter.before.burn, sep =";"))
      
      # Log current selected variables
      #cat("Current Selected Vars:", selected_vars, "\n")
      
      
      if(length(selected_vars) == 1){
        #for models with only one variable consider only additive nbd
        remaining_vars <- get_neighborhood(selected_vars, X,Y)
        
        # Additive neighborhood
        model_a_store <- NULL
        for (s in remaining_vars){
          p_addition <-  log_post_comparison(selected_vars,s)   # Compute the posterior probability
          model_a_store <- rbind(model_a_store, data.frame(p_addition, z = NA))
          model_a_store$z[nrow(model_a_store)] <- list(c(selected_vars,s))
        }
        
        total_posta <- logsumexp(model_a_store$p_addition)
        model_a_store$p_addition <- exp(model_a_store$p_addition - total_posta)
        accepted_models$add[nrow(accepted_models)] <- total_posta
        ad <- model_a_store[sample(nrow(model_a_store), 1, prob = (model_a_store$p_addition)^(1/temp[l])), ]#k+
        
        proposed_model <- ad #the additive model is proposed
        selected_vars <- unlist(proposed_model$z)
        accepted_models <- rbind(accepted_models, data.frame(post = proposed_model$p_addition,
                                                             model = NA,
                                                             add = NA,
                                                             rem = NA))
        accepted_models$model[nrow(accepted_models)] <- list(selected_vars)
        
        if (iter < foo){
          iter.before.burn <- iter.before.burn + 1
        }
        
        
        
      }else if(length(selected_vars) > 1){
        
        #for models having more than one variable
        remaining_vars <- get_neighborhood(selected_vars, X,Y)
        
        # Additive neighborhood and their posteriors
        model_a_store <- NULL
        for (s in remaining_vars){
          p_addition <-  log_post_comparison(selected_vars, s )   # Compute the posterior probability
          model_a_store <- rbind(model_a_store, data.frame(p_addition, z = NA))
          model_a_store$z[nrow(model_a_store)] <- list(c(selected_vars, s ))
        }
        total_posta <- logsumexp(model_a_store$p_addition)
        model_a_store$p_addition <- exp(model_a_store$p_addition - total_posta)
        accepted_models$add[nrow(accepted_models)] <- total_posta
        
        
        model_r_store <- NULL
        for (s in selected_vars){
          reduced_vars <- setdiff(selected_vars, s)
          #p_reduced <- post_comparison_re(selected_vars, s )   # Compute the posterior probability
          p_reduced <- -log_post_comparison(reduced_vars, s ) 
          model_r_store <- rbind(model_r_store, data.frame(p_reduced, z = NA))
          model_r_store$z[nrow(model_r_store)] <- list(reduced_vars)
        }
        
        total_postr <- logsumexp(model_r_store$p_reduced)
        model_r_store$p_reduced <- exp(model_r_store$p_reduced - total_postr)
        accepted_models$rem[nrow(accepted_models)] <- total_postr
        
        
        #sampling from additive and removed models
        ad <- model_a_store[sample(nrow(model_a_store), 1, 
                                   prob = (model_a_store$p_addition)^(1/temp[l])), ]#k+
        re <- model_r_store[sample(nrow(model_r_store), 1, 
                                   prob = (model_r_store$p_reduced)^(1/temp[l])), ]#k-
        
        
        add_prob <- (ad[,1])^(1/temp[l])
        rem_prob <- (re[,1])^(1/temp[l])
        prob_vec <- c(add_prob,rem_prob)
        #prob_vec <- prob_vec / sum(prob_vec)
        
        ##########################################
        #sampling from additive and remove model
        if(sample(1:2,1, prob = prob_vec) == 1){
          
          proposed_model <- ad #the additive model is proposed
          # accept-reject
          #remove neighborhood posterior for this proposed model
          total_postr_new <- 0
          for (s in unlist(proposed_model$z)){
            reduced_vars <- setdiff(unlist(proposed_model$z), s)
            total_postr_new <- total_postr_new -log_post_comparison(reduced_vars,s)
          }
          logR <- -proposed_model$p_addition - total_postr_new
          
          log_alpha <- min(0,logR)
          #ratio for accept-reject
          if(log(runif(1)) < log_alpha){
            selected_vars <- unlist(proposed_model$z)
            accepted_models <- rbind(accepted_models, 
                                     data.frame(post = proposed_model$p_addition,
                                                model = NA,
                                                add = NA,
                                                rem = NA))
            accepted_models$model[nrow(accepted_models)] <- list(selected_vars)
            if (iter < foo){
              iter.before.burn <- iter.before.burn + 1
            }
            
          }
        }else {
          
          proposed_model <- re#the removed model is proposed
          # accept-reject
          # additive neighborhood posterior for this proposed model
          remaining_vars <- get_neighborhood(unlist(proposed_model$z), X,Y)
          
          total_posta_new <-0
          for (s in remaining_vars){
            total_posta_new <- total_posta_new + log_post_comparison(unlist(proposed_model$z), s)
          }
          logR <- -proposed_model$p_reduced - total_posta_new
          
          
          log_alpha <- min(0, logR)
          if (log(runif(1)) < log_alpha) {
            selected_vars <- unlist(proposed_model$z)
            accepted_models <- rbind(accepted_models, data.frame(post = proposed_model$p_reduced,
                                                                 model = NA,
                                                                 add = NA,
                                                                 rem = NA))
            accepted_models$model[nrow(accepted_models)] <- list(selected_vars)
            if (iter < foo){
              iter.before.burn <- iter.before.burn + 1
            }
          }
          
        }
        
        # Stop if the number of selected variables reaches p-1
        if (length(selected_vars) == p - 1) {
          break
        }
      }
    }#inner loop ends
    #print(accepted_models[,1:2])###for checking all outputs
    #print(proposed_model$z)# if you want to check each acceptance
  }#outer loop ends
  
  toc <- proc.time()[3]
  a <- toc-tic
  print(a)
  
  
  #burn-in
  answer <- accepted_models[-(1:iter.before.burn),1:2]
  
  #MAP model
  map_indices <- which(answer$post == max(answer$post))
  map_model <- sort(answer$model[[map_indices[which.min(sapply(answer$model[map_indices], length))]]])
  
  #Mode model
  b <- sort(table(sapply(answer$model, function(x) toString(sort(x)))), decreasing = TRUE)
  s <- names(b[1])
  mode_model <- as.integer(unlist(strsplit(s, ",")))
  
  
  #Median model
  prob <- table(unlist(answer$model))/dim(answer)[1]
  median_model <- as.numeric(names(prob[prob>0.5]))
  
  #Master predictor
  star_matrix <- matrix(0, nrow = p, ncol = q)
  for(index in 1:dim(answer)[1]) {
    current_model <- answer$model[[index]] 
    model_post <- answer$post[[index]]  
    B_gamma_mat <- B_gamma(current_model)
    
    for(i in seq_along(current_model)) {
      j <- current_model[i]
      star_matrix[j, ] <- star_matrix[j, ] + model_post * B_gamma_mat[i, ]
    }
  }
  
  # Compute alpha_j for each predictor j:
  MPI <- apply(star_matrix, 1, function(row_vals) {
    prod_val <- prod(abs(row_vals))  # Multiply absolute values over l=1 to q
    tanh(prod_val/2)
  })
  
  return(kutu = list(
    MAP = map_model,
    Mode = mode_model,
    MPM = median_model,
    Star = data.frame(
    predictor = which(MPI != 0),
    MPI = MPI[MPI != 0]
    )))
}


####################################################
#competitors


MCC <- function(X, Y, threshold = NULL) {
  
  # Estimate sample covariance matrix for Y
  cov_Y <- cov(Y)
  inv_cov_Y <- solve(cov_Y)
  
  # Compute multiple correlation coefficient for each predictor
  mcc_values <- sapply(1:p, function(k) {
    xk <- X[, k]  # Extract the k-th predictor
    gamma_k <- inv_cov_Y %*% colMeans(Y * xk)  # Compute gamma_k (q x 1 vector)
    rho_k2 <- mean((Y %*% gamma_k) * xk)  # Compute rho_k^2
    return(rho_k2)
  })
  
  # Threshold for screening
  dn <- ceiling(n / log(n))
  sorted_indices <- order(mcc_values, decreasing = TRUE)
  selected_predictors <- sorted_indices[1:dn]
  return(selected_predictors)
}  # Custom MCC


##################################################
results_list <- list()
# g-values 
x_values <- c(1,2,2.5)

for (setup_id in 1:4) {
  s <- cov_setups[[setup_id]]
  i = setup_id
  s <- cov_setups[[setup_id]]
  n <- s$n
  p <- s$p
  p0 <- s$p0
  q <- s$q
  
  
  cat("Running setup", i, "\n")

  
  prefix <- paste0( p0, "_", q, "_", setup_id)
  # Load data
  X <- as.matrix(read.csv(paste0("X_", prefix, ".csv"), header = TRUE))
  Y <- as.matrix(read.csv(paste0("Y_", prefix, ".csv"), header = TRUE))
  B <- as.matrix(read.csv(paste0("B_", prefix, ".csv"), header = TRUE))
  E <- Y - X %*% B
  indx.beta <- as.vector(read.csv(paste0("S_", prefix, ".csv"), header = TRUE))$x

  x_values <- c(1,2,2.5)
  g_values <- round(c(log(p^x_values),log(q^x_values), log(p)*log(q), exp(n^(x_values/10)),p*(1:10)/10,q^(q/2),p^2))
  
  
  snr <- sum(diag(t(B[indx.beta,]) %*%cov(X[,indx.beta])%*% B[indx.beta,]))/sum(diag(cov(E)))
  
  # --- Elastic Net (glmnet) ---
  fit_enet <- lapply(1:q, function(j) {
    cv.glmnet(X, Y[, j], alpha = 0.5)
  })
  enet_mat <- sapply(fit_enet, function(f) coef(f, s = "lambda.min")[-1])
  row_norms <- apply(enet_mat, 1, function(r) sqrt(sum(r^2)))
  enet <- order(row_norms, decreasing = TRUE)[1:ceiling(n/log(n))]
  
  # --- SCAD (ncvreg) ---
  fit_scad <- lapply(1:q, function(j) {
    cv.ncvreg(X, Y[, j], penalty = "SCAD")
  })
  scad_mat <- sapply(fit_scad, function(f) coef(f, s = "lambda.min")[-1])
  row_norms <- apply(scad_mat, 1, function(r) sqrt(sum(r^2)))
  scad <- order(row_norms, decreasing = TRUE)[1:ceiling(n/log(n))]
  
  
  # --- BLASSO --
  #fit_BayesSUR <- BayesSUR(
  #  Y = Y, X = X, outFilePath = "BayesSUR_output", 
  #  covariancePrior = "IW",
  #  gammaPrior = "hierarchical",
  #  nIter = 1000, 
  #  burnin = 500, 
  #  nChains = 1
  #)
  #pi_vec <- read.table("BayesSUR_output/data_dSUR_pi_out.txt", header = FALSE)[[1]]
  #bayessur <- order(pi_vec, decreasing = TRUE)[1:ceiling(n/log(n))]
  
  # --- MCC ---
  mcc <- MCC(X,Y,NULL)
  
  # --- MBSP ---
  mbsp_model <- MBSP(Y=Y, X=X, max_steps=1500, burnin=1000, model_criteria=FALSE)
  
  # our method
  results <- data.frame()
  STAR <- list()
  for (g in g_values) {
    result <- s_53_log(X, Y = Y, L = 20, J = 75, burnin = 0.6666,
                       temp = seq(1, 0.1, length.out = 100), g = g)
    
    res_map <- cbind(g = g, Method = "MAP", Predictors = paste(result$MAP, collapse = ","),
                     as.data.frame(PRFJ(result$MAP, indx.beta)))
                     
    res_mpm <- cbind(g = g, Method = "MPM", Predictors = paste(result$MPM, collapse = ","),
                     as.data.frame(PRFJ(result$MPM, indx.beta)))
    res_mode <- cbind(g = g, Method = "Mode", Predictors = paste(result$Mode, collapse = ","),
                      as.data.frame(PRFJ(result$Mode, indx.beta)))
    
    results <- rbind(results, res_map, res_mpm, res_mode)
    
    
    STAR[[as.character(g)]] <- result$Star
  }
  
  results_list[[i]] <- list(
    cov_y_type = s$cov_y,
    cov_x_type = s$cov_x,
    snr = snr,
    real_B = B,
    true_predictor = indx.beta,
    result = results,
    star = STAR,
    comp_mbsp = cbind(Predictors = paste(mbsp_model$active_predictors, collapse = ","),
                      as.data.frame(PRFJ(mbsp_model$active_predictors, indx.beta))),
    comp_bsrrr = NA,#will be updated separately
    #cbind(Predictors = paste(comp_bsrrr,collapse = ","),as.data.frame(PRFJ(comp_bsrrr, indx.beta)),
    
    comp_mcc = cbind(Predictors = paste(mcc, collapse = ","),
                     as.data.frame(PRFJ(mcc, indx.beta))),
    
    comp_bayessur = NA,#cbind(Predictors = paste(bayessur, collapse = ","),
    #                      as.data.frame(PRFJ(bayessur, indx.beta))),
    
    comp_scad = cbind(Predictors = paste(scad, collapse = ","),
                      as.data.frame(PRFJ(scad, indx.beta))),
    comp_enet = cbind(Predictors = paste(enet, collapse = ","),
                      as.data.frame(PRFJ(enet, indx.beta)))
                      
    
  )
  
}

###################################################
#parallel
library(foreach)
library(doParallel)

cores <- parallel::detectCores() - 2  # Leave one core free
cl <- makeCluster(cores)
registerDoParallel(cl)

setup_grid <- expand.grid(
  setup_id = 1:4
)


results_list <- foreach(setup_id = 1:4, .packages = c("glmnet", "ncvreg", "MASS", "monomvn","MBSP",
                                                      "BayesSUR","LaplacesDemon","matrixsampling")) %dopar% {
  s <- cov_setups[[setup_id]]
  i <- setup_id
  n <- s$n
  p <- s$p
  p0 <- s$p0
  q <- s$q
  
  cat("Running setup", i, "\n")
  
  prefix <- paste0(p0, "_", q, "_", setup_id)
  X <- as.matrix(read.csv(paste0("X_", prefix, ".csv"), header = TRUE))
  Y <- as.matrix(read.csv(paste0("Y_", prefix, ".csv"), header = TRUE))
  B <- as.matrix(read.csv(paste0("B_", prefix, ".csv"), header = TRUE))
  E <- Y - X %*% B
  indx.beta <- as.vector(read.csv(paste0("S_", prefix, ".csv"), header = TRUE))$x
  
  
  g_values <- round(c(log(p^x_values), log(q^x_values), log(p)*log(q), exp(n^(x_values/10)), p*(1:10)/10, q^(q/2), p^2))
  snr <- sum(diag(t(B[indx.beta,]) %*% cov(X[, indx.beta]) %*% B[indx.beta,])) / sum(diag(cov(E)))
  
  # Elastic Net
  fit_enet <- lapply(1:q, function(j) {
    cv.glmnet(X, Y[, j], alpha = 0.5)
  })
  enet_mat <- sapply(fit_enet, function(f) coef(f, s = "lambda.min")[-1])
  row_norms <- apply(enet_mat, 1, function(r) sqrt(sum(r^2)))
  enet <- order(row_norms, decreasing = TRUE)[1:ceiling(n/log(n))]
  
  # SCAD
  fit_scad <- lapply(1:q, function(j) {
    cv.ncvreg(X, Y[, j], penalty = "SCAD")
  })
  scad_mat <- sapply(fit_scad, function(f) coef(f, s = "lambda.min")[-1])
  row_norms <- apply(scad_mat, 1, function(r) sqrt(sum(r^2)))
  scad <- order(row_norms, decreasing = TRUE)[1:ceiling(n/log(n))]
  
  # MCC
  mcc <- MCC(X, Y, NULL)
  
  # MBSP
  mbsp_model <- MBSP(Y = Y, X = X, max_steps = 1500, burnin = 1000, model_criteria = FALSE)
  
  # Our method
  results <- data.frame()
  STAR <- list()
  for (g in g_values) {
    result <- s_53_log(X, Y = Y, L = 20, J = 75, burnin = 0.6666,
                       temp = seq(1, 0.1, length.out = 100), g = g)
    
    res_map <- cbind(g = g, Method = "MAP", Predictors = paste(result$MAP, collapse = ","),
                     as.data.frame(PRFJ(result$MAP, indx.beta)))
    res_mpm <- cbind(g = g, Method = "MPM", Predictors = paste(result$MPM, collapse = ","),
                     as.data.frame(PRFJ(result$MPM, indx.beta)))
    res_mode <- cbind(g = g, Method = "Mode", Predictors = paste(result$Mode, collapse = ","),
                      as.data.frame(PRFJ(result$Mode, indx.beta)))
    
    results <- rbind(results, res_map, res_mpm, res_mode)
    STAR[[as.character(g)]] <- result$Star
  }
  
  list(
    cov_y_type = s$cov_y,
    cov_x_type = s$cov_x,
    snr = snr,
    real_B = B,
    true_predictor = indx.beta,
    result = results,
    star = STAR,
    comp_mbsp = cbind(Predictors = paste(mbsp_model$active_predictors, collapse = ","),
                      as.data.frame(PRFJ(mbsp_model$active_predictors, indx.beta))),
    comp_bsrrr = NA,
    comp_mcc = cbind(Predictors = paste(mcc, collapse = ","),
                     as.data.frame(PRFJ(mcc, indx.beta))),
    comp_scad = cbind(Predictors = paste(scad, collapse = ","),
                      as.data.frame(PRFJ(scad, indx.beta))),
    comp_enet = cbind(Predictors = paste(enet, collapse = ","),
                      as.data.frame(PRFJ(enet, indx.beta)))
  )
}

names(results_list) <- paste0("setup_", 1:4)

save(results_list, file = "simulation_results_all.RData")

stopCluster(cl)
library(xtable)
data <- results_list[[1]]$result[-3]
latex <- xtable(data)

for (i in 1:22){
  results_list[["setup_1"]][["star"]][[i]] <- 
    cbind(results_list[["setup_1"]][["star"]][[i]],(results_list[["setup_1"]][["star"]][[i]]$MPI +1)/2)
}


################################################################################
#BSRRR
library(grpreg)

# Load C_mean p × q)
indx.betas <- list()
for (i in 1:6){
  C_mean <- as.matrix(read.csv(paste("C_bsrrr_",i,".csv", sep=""), header = FALSE))
  p <- nrow(C_mean)
  q <- ncol(C_mean)
  
  # Flatten C_mean column-wise ??? response vector
  y_vec <- as.vector(C_mean)
  
  # Create design matrix Z: p*q × p
  Z <- matrix(0, nrow = p*q, ncol = p)
  for (j in 1:q) {
    rows <- ((j - 1) * p + 1):(j * p)
    Z[rows, ] <- diag(p)
  }
  
  # Group index: each column in Z is a group (row of C_mean)
  group <- 1:p
  
  # Cross-validated fit
  cvfit <- cv.grpreg(Z, y_vec, group = group, penalty = "grLasso")
  
  # Extract selected groups at lambda.min
  coef_vec <- coef(cvfit, lambda = cvfit$lambda.min)[-1]  # drop intercept
  group_norms <- sapply(split(coef_vec, group), function(g) sqrt(sum(g^2)))
  bsrrr <- order(group_norms, decreasing = TRUE)[1:(n/log(n))]
  
  B <- results_list[[i]][["real_B"]]
  indx.beta <- unique(which(B!=0, arr.ind = TRUE)[,1])
  results_list[[i]][["comp_bsrrr"]] = cbind(Predictors = paste(bsrrr, collapse = ","),
                                            as.data.frame(PRFJ(bsrrr, indx.beta)),
                                            Bselected_B0=norm(B[bsrrr,],"2")/norm(B[indx.beta,]))
  
}

