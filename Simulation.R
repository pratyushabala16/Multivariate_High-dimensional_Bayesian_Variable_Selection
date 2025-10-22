rm(list=ls())
set.seed(123)  # for reproducibility
library(Matrix)
library(MASS)
library(LaplacesDemon) # sampling omega
library(matrixsampling)# or LaplacesDemon to sample MxVT
###############################################
# hyperparameter tau
tau <- 1e-5
##################################################################
#Our functions
### log posterior function
log_post <- function(X,Y,selected_vars){
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
log_post_comparison <- function(X,Y,selected_vars, s){
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
  
### sampling for posterior mean of B from MxVT(M,Sigma,Omega,neu)
B_gamma <- function(X,Y,selected_vars){
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

# Main function
multivariate_hd_bayesian_vs <- function(X, Y, temp = seq(1, 0.1, length.out = 20) , 
                     L = 20,J=20, burnin = 0.5,g = g, M =n/log(n)){
  
  #scaling
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  #tau <- 1e-5
  
  
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
    p_addition <- log_post(X,Y,s)  # Compute the marginal posterior probability
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
      
      # Log current selected variables if you want
      #cat("Current Selected Vars:", selected_vars, "\n")
      
      
      if(length(selected_vars) == 1){
        #for models with only one variable consider only additive nbd
        remaining_vars <- get_neighborhood(selected_vars, X,Y)
        
        # Additive neighborhood
        model_a_store <- NULL
        for (s in remaining_vars){
          p_addition <-  log_post_comparison(X,Y,selected_vars,s)   # Compute the log posterior probability
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
          p_addition <-  log_post_comparison(X,Y,selected_vars, s )   # Compute the posterior probability
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
          p_reduced <- -log_post_comparison(X,Y,reduced_vars, s ) 
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
        
        ##########################################
        #sampling from additive and remove model
        if(sample(1:2,1, prob = prob_vec) == 1){
          
          proposed_model <- ad #the additive model is proposed
          # accept-reject
          # remove neighborhood posterior for this proposed model
          total_postr_new <- 0
          for (s in unlist(proposed_model$z)){
            reduced_vars <- setdiff(unlist(proposed_model$z), s)
            total_postr_new <- total_postr_new -log_post_comparison(X,Y,reduced_vars,s)
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
            total_posta_new <- total_posta_new + log_post_comparison(X,Y,unlist(proposed_model$z), s)
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
  for(index in 1:nrow(answer)) {
    current_model <- answer$model[[index]] 
    model_post <- answer$post[[index]]  
    B_gamma_mat <- B_gamma(X,Y,current_model)
    
    for(i in seq_along(current_model)) {
      j <- current_model[i]
      star_matrix[j, ] <- star_matrix[j, ] + (model_post) * B_gamma_mat[i, ]
    }
  }
  star_matrix <- star_matrix/nrow(answer)
  
  # Compute alpha_j for each predictor j:
  MPI <- apply(star_matrix, 1, function(row_vals) {
    prod_val <- prod(abs(row_vals))  # Multiply absolute values over l=1 to q
    a <-(prod_val^(1/q))
    return(tanh(a/2))# psi function
  })
  predictors <- which(MPI != 0)
  MPI <- MPI[predictors]
  MPI <- round((MPI-min(MPI))/(max(MPI)-min(MPI)),3)
  Star = data.frame(predictors,MPI)
  return(kutu = list(
    MAP = map_model,
    Mode = mode_model,
    MPM = median_model,
    Star =Star
    ))
}


### post processing functions for comparison
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

##################################################
# Simulation Run
# data
n <- 200 # no of samples
p <- 1000 # no of covariates
p0 <- 15 # no of true covariates
q <- 10 # no of responses; q < p0
                         
X = matrix(rnorm(n*p),n,p)
X = scale(X) # cov_x = "indpt"

# simulate B covariate matrix of your own choice
B <- matrix(0, nrow = p, ncol = q) 
true_indices <- sort(sample(1:p, p0)) # all imporatant covariates
top <- strong_rows[1:floor(p0/3)] # highly important covariates
mid <- strong_rows[ceiling(p0/3):floor(2*p0/3)] # mid important covariates
low <- strong_rows[ceiling(2*p0/3):p0] # low important covariates

B[top, ] <- matrix(sample(c(-2, -1, 1, 2), length(top)*q, prob=c(0.1,0.4,0.4,0.1),replace = TRUE), nrow = length(top))
for (r in mid) {
    pos <- sample(1:q, q/2)
    B[r, pos] <- sample(c(-1,-2, 2,1),prob=c(0.2,0.3,0.3,0.2), length(pos), replace = TRUE)
}
for (r in low) {
    pos <- sample(1:q, ceiling(q/4))
    B[r, pos] <- sample(c(-2, 2), length(pos), replace = TRUE)
}
# adjust for SNR
norms_q <- NULL
for ( i in strong_rows){
    norms_q <- rbind(norms_q, norm(B[i,1:q],type = '2'))
}
B <- B/floor(sqrt(max(norms_q)))

# simualate y
rho =0.7 # if y has dense covariance 
cov_y = (rho * matrix(1, q, q) + (1 - rho) * diag(q))
E <- mvrnorm(n, rep(0, q), cov_y)
Y <- scale(X %*% B + E)                        
# hyperparameter g-values 
x_values <- c(1,2,2.5)  
g_values <- round(c(log(p^x_values),log(q^x_values))) # should try other values too
names(g_values) <- c("log(p)","2log(p)","2.5log(p)","log(q)","2log(q)","2.5log(q)")
results <- data.frame() # stores the accepted models and metrics
STAR <- list() # stores the MPP values
for (g in g_values) {
    result <-multivariate_hd_bayesian_vs(X, Y = Y, L = 50, J = 100, burnin = 0.4,
                       temp = seq(1, 0.1, length.out = 100), g = g)
    
    res_map <- cbind(g = g, Method = "MAP", Predictors = paste(result$MAP, collapse = ","),
                     as.data.frame(PRFJ(result$MAP, indx.beta)))
                     
    res_mpm <- cbind(g = g, Method = "MPM", Predictors = paste(result$MPM, collapse = ","),
                     as.data.frame(PRFJ(result$MPM, indx.beta)))
    res_mode <- cbind(g = g, Method = "Mode", Predictors = paste(result$Mode, collapse = ","),
                      as.data.frame(PRFJ(result$Mode, indx.beta)))
    
    results <- rbind(results, res_map, res_mpm, res_mode)
    STAR[[as.character(g)]] <- result$Star
    names(STAR) <- names(g_values) 
}
save(results,file="results.Rdata")


