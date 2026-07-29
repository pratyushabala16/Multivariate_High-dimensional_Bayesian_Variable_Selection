rm(list=ls())
set.seed(123)  # for reproducibility
library(Matrix)
library(MASS)
library(LaplacesDemon) # sampling omega
library(matrixsampling)# or LaplacesDemon to sample MxVT

##################################################################
# posterior mean for B
B_gamma <- function(X,Y,selected_vars){
  ### from MxVT(M,Sigma,Omega,neu)
  a <- 1/(1+ 1/g)
  b <- 1+ 1/tau
  c <- q + 1
  pgamma <- length(selected_vars)
  v_n <- 1-b^{-1}
  u_n <- a^{-1} - b^{-1}
  x <- X[, selected_vars]
  P <- x%*%solve(t(x)%*%x)%*%t(x)
  # posterior mean for B_gamma
  U <- solve(t(x)%*%x)/u_n
  M <- v_n*U%*%t(x)%*%Y
  #V <- diag(q)
  #V <- {v_n*t(Y)%*%(diag(n)-(v_n/u_n)*P)%*%Y +V}/(n+c-q+1)
  #V_fixed <- as.matrix(nearPD(V)$mat)
  #library(matrixsampling)
  #B_g <- rmatrixt(n=10, nu=(n+c-q+1), M=M, U=U, V=V_fixed, checkSymmetry = FALSE, keep = TRUE)
  #B <-  apply(B_g, c(1, 2), mean)
  return(M)
}

# --------------- MAIN FUNCTION -------------------------------
### multivariate Bayesian high-dimensional variable selection
mbhvs <- function(X, Y, L = 20, J=20, tau = tau,
                  temp = seq(1, 0.1, length.out = 20),#L , 
                  burnin = 0.4, g = g, M = n/log(n)){
  
  # scaling
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  
  # Initialize empty set of selected variables
  selected_vars <- integer(0)
  model_a_store <- NULL
  model_r_store <- NULL
  
  # neighborhood function
  get_neighborhood <- function(selected_vars, X, Y) {
    Z <- matrix(1, nrow(X), 1)
    
    if (length(selected_vars) > 0) {
      Z <- cbind(Z, X[, selected_vars, drop = FALSE])
    }
    
    Rmat <- qr.resid(qr(Z), Y)
    
    h <- crossprod(Rmat, X)
    h_norms <- sqrt(colSums(h^2))
    
    top_vars <- order(h_norms, decreasing = TRUE)[
      seq_len(min(M, ncol(X)))
    ]
    
    setdiff(top_vars, selected_vars)
  }
  
  logsumexp <- function(log_vals) {
    m = max(log_vals)
    return(m + log(sum(exp(log_vals - m))))
  }
  

  #Rcpp::sourceCpp("mbhvs_cpp.cpp")
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
    p_addition <- log_post_cpp(X,Y,as.integer(s), g, tau)  # Compute the marginal posterior probability
    model_a_store <- rbind(model_a_store, data.frame(p_addition,s))
  }
  selected_vars <- model_a_store$s[which.max(model_a_store$p_addition)]#map
  accepted_models <- data.frame(post = max(model_a_store$p_addition),
                                model = selected_vars,
                                add=NA,
                                rem=NA,
                                iter=0)
  
  
  iter <- 0
  iter.before.burn <- 0
  
  foo <- ceiling((L*J)*burnin)
  # Stores the model occupied at every iteration
  #visited_models <- vector("list", L * J)
  #visited_logpost <- rep(NA, L * J)
  # Optional: records whether that iteration accepted a move
  move_accepted <- rep(NA, L * J)
  
  
  #starting the loop
  for (l in 1:L) {
    for (i in 1:J){
      iter <- (l - 1) * J + i
      #print(paste(l,i,iter,iter.before.burn, sep =";"))
      if (iter %% 100 == 0) {
        message("Iteration: ", iter)
      }
      accepted_this_iter <- FALSE
      # Log current selected variables
      #cat("Current Selected Vars:", selected_vars, "\n")
      
      
      if(length(selected_vars) == 1){
        #for models with only one variable consider only additive nbd
        remaining_vars <- get_neighborhood(selected_vars, X,Y)
        
        # Additive neighborhood
        model_a_store <- NULL
        for (s in remaining_vars){
          p_addition <-  log_post_comparison_cpp(
            X, Y,
            as.integer(selected_vars),
            as.integer(s),
            g, tau
          )   # Compute the log posterior probability
          model_a_store <- rbind(model_a_store, data.frame(p_addition, z = NA))
          model_a_store$z[nrow(model_a_store)] <- list(c(selected_vars,s))
        }
        
        total_posta <- logsumexp(model_a_store$p_addition)
        model_a_store$p_addition <- exp(model_a_store$p_addition - total_posta)
        accepted_models$add[nrow(accepted_models)] <- total_posta
        ad <- model_a_store[sample(nrow(model_a_store), 1, prob = (model_a_store$p_addition)^(1/temp[l])), ]#k+
        
        proposed_model <- ad #the additive model is proposed
        selected_vars <- unlist(proposed_model$z)
        accepted_this_iter <- TRUE
        accepted_models <- rbind(accepted_models, data.frame(post = proposed_model$p_addition,
                                                             model = NA,
                                                             iter = iter,
                                                             add = NA,
                                                             rem = NA))
        accepted_models$model[nrow(accepted_models)] <- list(selected_vars)
        #visited_models[[iter]] <- selected_vars
        #visited_logpost[iter] <- accepted_models$post[accepted_models$iter == iter]
        
        
        if (iter < foo){
          iter.before.burn <- iter.before.burn + 1
        }
        
        
        
      }else if(length(selected_vars) > 1){
        
        #for models having more than one variable
        remaining_vars <- get_neighborhood(selected_vars, X,Y)
        
        # Additive neighborhood and their posteriors
        model_a_store <- NULL
        for (s in remaining_vars){
          p_addition <-  log_post_comparison_cpp(
            X, Y,
            as.integer(selected_vars),
            as.integer(s),
            g, tau
          )  # Compute the posterior probability
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
          p_reduced <- -log_post_comparison_cpp(
            X, Y,
            as.integer(reduced_vars),
            as.integer(s),
            g, tau
          )
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
            total_postr_new <- total_postr_new -log_post_comparison_cpp(
              X, Y,
              as.integer(reduced_vars),
              as.integer(s),
              g, tau
            )
          }
          logR <- -proposed_model$p_addition - total_postr_new
          
          log_alpha <- min(0,logR)
          #ratio for accept-reject
          if(log(runif(1)) < log_alpha){
            selected_vars <- unlist(proposed_model$z)
            accepted_this_iter <- TRUE
            accepted_models <- rbind(accepted_models, 
                                     data.frame(post = proposed_model$p_addition,
                                                model = NA,
                                                add = NA,
                                                rem = NA,
                                                iter = iter))
            accepted_models$model[nrow(accepted_models)] <- list(selected_vars)
            if (iter < foo){
              iter.before.burn <- iter.before.burn + 1
            }
            
          }
        }else{
          
          proposed_model <- re#the removed model is proposed
          # accept-reject
          # additive neighborhood posterior for this proposed model
          remaining_vars <- get_neighborhood(unlist(proposed_model$z), X,Y)
          
          total_posta_new <-0
          for (s in remaining_vars){
            total_posta_new <- total_posta_new + log_post_comparison_cpp(
              X,Y, as.integer(unlist(proposed_model$z)), as.integer(s), g, tau)
          }
          logR <- -proposed_model$p_reduced - total_posta_new
          
          
          log_alpha <- min(0, logR)
          if (log(runif(1)) < log_alpha) {
            selected_vars <- unlist(proposed_model$z)
            accepted_this_iter <- TRUE
            accepted_models <- rbind(accepted_models, data.frame(post = proposed_model$p_reduced,
                                                                 model = NA,
                                                                 add = NA,
                                                                 rem = NA,
                                                                 iter = iter))
            accepted_models$model[nrow(accepted_models)] <- list(selected_vars)
            if (iter < foo){
              iter.before.burn <- iter.before.burn + 1
            }
          }
          
        }
        # Save the state occupied at this iteration.
        #visited_models[[iter]] <- selected_vars
        #visited_logpost[iter] <- accepted_models$post[nrow(accepted_models)]
        move_accepted[iter] <- accepted_this_iter
        # Stop if the number of selected variables reaches p-1
        # max_model_size <- min(p, n)
        if (length(selected_vars) == p - 1) {
          break
        }
      }
    

      #visited_logpost[iter] <- accepted_models$post[accepted_models$iter == iter]

      }#inner loop ends
    
    #print(accepted_models[,1:2])###for checking all outputs
    #print(proposed_model$z)# if you want to check each acceptance
  }#outer loop ends
  
  toc <- proc.time()[3]
  time_needed <- toc-tic
  print(time_needed)
  
  # acceptance rate
  ac_rate <- mean(move_accepted[(foo + 1):(L * J)], na.rm = TRUE)

  #post_models <- visited_models[(foo + 1):(L * J)]
  #post_models <- post_models[!vapply(post_models, is.null, logical(1))]
  #post_logpost <- visited_logpost[(foo + 1):(L * J)]
  
  #MAP model
  #map_indices <- which(post_logpost == max(post_logpost))
  ## length of models with highest posterior
  #len <- sapply(post_models, length)[map_indices]
  #min_indx <- which.min(len)
  #Map_model <- sort(post_models[map_indices][[min_indx]])
  
  
  # Mode model: most frequently visited exact model
  #model_key <- vapply(
  #  post_models,
  #  function(z) paste(sort(z), collapse = ","),
  #  character(1)
  #)
  
  #mode_key <- names(which.max(table(model_key)))
  #Mode_model <- as.integer(strsplit(mode_key, ",", fixed = TRUE)[[1]])
  
  
  # Median probability model: variables with inclusion frequency >= 0.5
  #PIP <- sapply(seq_len(p), function(j) {
  #  mean(vapply(post_models, function(z) j %in% z, logical(1)))
  # })
  
  #Median_model <- which(PIP >= 0.5)
  
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
    current_model <- answer$model[[index]] #post_models[[index]]
    model_post <- answer$post[[index]] #post_logpost[index]
    B_gamma_mat <- B_gamma(X,Y,current_model)
    
    for(i in seq_along(current_model)) {
      j <- current_model[i]
      star_matrix[j, ] <- star_matrix[j, ] + B_gamma_mat[i, ]*(model_post) 
    }
  }
  
  star_matrix <- star_matrix/nrow(answer)
  
  # Compute alpha_j for each predictor j:
  MPI <- apply(star_matrix, 1, function(row_vals) {
    prod_val <- prod(abs(row_vals))  # Multiply absolute values over l=1 to q
    a <-(prod_val^(1/q))
    return(tanh(a/2))#psi
  })
  

  predictors <- which(MPI != 0)
  MPI <- round(MPI[predictors],3)
  MPI <- round((MPI-min(MPI))/(max(MPI)-min(MPI)),3)
  Star = data.frame(predictors,MPI)
  return(kutu = list(
    models = answer$model,
    logposterior = answer$post,
    MAP = map_model,
    Mode = mode_model,
    MPM = median_model,
    Star = Star,
    time_needed = time_needed,
    ac_rate = ac_rate
  ))
}


# ------------------Before using Main function ----------------
Rcpp::sourceCpp("mbhvs_cpp.cpp")

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
  
  return(list(precision = round(precision,2),
              recall = round(recall,2),
              F1 = round(F1,2),
              Jaccard_index = round(Jaccard,2),
              norm_BbyB0 = round(norm_BbyB0,2)
  ))
}
