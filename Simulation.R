###FINAL CODE

#      
#simulation generate true model
rm(list = ls())
set.seed(123)
p =1000
q =3
n= 100
# for |gamma_0| = 5
indx.beta = c(1,4,33,78,99)
xd0 = rep(0,p);xd0[indx.beta]=1
B = matrix(0, nrow = p, ncol = q)
bt0 = rep(0,p);
bt0[indx.beta]=c(1,1.25,1.5,1.75,2)*sample(c(1,-1),5,replace=TRUE)
B[,1]=bt0
bt0 = rep(0,p);
bt0[c(1,4,33,78,99,9)]=c(2,1.75,1.25,1,1.75,1.25)*sample(c(1,-1),6,replace=TRUE)
B[,2]=bt0
bt0 = rep(0,p);
bt0[c(1,4,33,78,99,5,9)]=c(1.75,1,1.3,1.5,1.25,2,1.3)*sample(c(1,-1),7,replace=TRUE)
B[,3]=bt0

##putting noise
bt0[c(5,9,120,490,659,200,400)]=c(1.75,1,1.3,1.5,1.25,2,1.3)*sample(c(1,-1),7,replace=TRUE)
B[,2]=bt0
all <- c(indx.beta,5,9,120,490,659,200,400)
B[all,]
B_small <- B


# for |gamma_0| = 10 put extra significant covariates
indx.beta2 = c(300,600,479,500,800)
B[indx.beta2,1]=c(1,1.25,1.5,1.75,2)*sample(c(1,-1),5,replace=TRUE)
B[indx.beta2,2]=c(2,1.75,1.25,1,1.75)*sample(c(1,-1),5,replace=TRUE)
B[indx.beta2,3]=c(1.75,1.3,1.5,1.25,2)*sample(c(1,-1),5,replace=TRUE)
B[c(all,indx.beta2),]

##Independent X
X = matrix(rnorm(n*p),n,p)
X = scale(X)


##AR-1 dependency 
times <- 1:p
rho <- 0.5
H <- abs(outer(times, times, "-"))
V <- rho^H
mu <- rep(0, p)
# Rows of X are simulated from MVN(0,V)
X <- mvtnorm::rmvnorm(n, mu, V)
X <- scale(X, center=TRUE, scale=TRUE)


#block dependence inside X
block_1 <- length(all)              # Block 1: True predictors and related
block_2 <- p - block_1              # Block 2: Remaining predictors

# Generate covariance for Block 1
rho_1 <- 0.7  # Correlation within Block 1
Sigma_block1 <- matrix(rho_1, nrow = block_1, ncol = block_1) + diag(1 - rho_1, block_1)

# Generate covariance for Block 2
rho_2 <- 0.4  # Correlation within Block 2
Sigma_block2 <- matrix(rho_2, nrow = block_2, ncol = block_2) + diag(1 - rho_2, block_2)

# Decide dependency between blocks (optional)
rho_between <- 0  # Correlation between blocks (set to 0 for independence)

# Combine blocks into a full covariance matrix
Sigma_X <- matrix(0, nrow = p, ncol = p)  # Initialize
Sigma_X[all, all] <- Sigma_block1
remaining_indices <- setdiff(1:p, all)
Sigma_X[remaining_indices, remaining_indices] <- Sigma_block2
Sigma_X[all, remaining_indices] <- rho_between
Sigma_X[remaining_indices, all] <- rho_between
# Check positive definiteness (adjust if needed)
if (!all(eigen(Sigma_X)$values > 0)) {
  diag(Sigma_X) <- diag(Sigma_X) + 1e-6  # Add small value to diagonal for stability
}

# Simulate X
mu <- rep(0, p)  # Mean vector
X <- mvtnorm::rmvnorm(n = n, mu,Sigma_X)
X <- scale(X, center = TRUE, scale = TRUE)



# generating responses
cros <- crossprod(t(X),B) 
y = cros + rnorm(n)*sqrt(1.5)
y = scale(y, center=TRUE, scale=TRUE)
Y<- y

# error covariance
A = matrix(0.15, nrow = 3, ncol = 3) + diag(1 - 0.15, 3)
y1 = cros + mvtnorm::rmvnorm(n,rep(0,3), sqrt(1.5)*A)
y1=scale(y1)
###check
norm(X[, indx.beta]%*%B[indx.beta,], type = '2') > norm(X[,-(indx.beta)]%*%B[-(indx.beta),], type = '2')
norm(X[, indx.beta]%*%B[indx.beta,], type = 'F') > norm(X[,-(indx.beta)]%*%B[-(indx.beta),], type = 'F')

constant <- norm(X[, indx.beta]%*%B[indx.beta,], type = '2') / norm(X[,-(indx.beta)]%*%B[-(indx.beta),], type = '2')

## constant2 for |gamma_0|=10
norm(X[, c(indx.beta,indx.beta2)]%*%B[c(indx.beta,indx.beta2),], type = '2') > norm(X[,-(c(indx.beta,indx.beta2))]%*%B[-(c(indx.beta,indx.beta2)),], type = '2')
norm(X[, c(indx.beta,indx.beta2)]%*%B[c(indx.beta,indx.beta2),], type = 'F') > norm(X[,-(c(indx.beta,indx.beta2))]%*%B[-(c(indx.beta,indx.beta2)),], type = 'F')

constant2 <- norm(X[, c(indx.beta,indx.beta2)]%*%B[c(indx.beta,indx.beta2),], type = '2') / norm(X[,-(c(indx.beta,indx.beta2))]%*%B[-(c(indx.beta,indx.beta2)),], type = '2')

# hyperparameters
tau <- 1e-5
g <- p^2/2##need to verify each time

# posterior probability
post<- function(selected_vars, c=3){
  pgamma <- length(selected_vars)
  V <- diag(q)
  
  #rho = 0.5
  #V = matrix(0, nrow = q, ncol = q)
  
  # Fill in AR(1) covariance
  #for (i in 1:q) {
   # for (j in 1:q) {
    #  V[i, j] = rho^abs(i - j)
    #}
  #}
  
  a <- 1/(1+ 1/g)
  b <- 1+ 1/tau
  x <- X[, selected_vars]
  P <- x%*%solve(t(x)%*%x)%*%t(x)
  zai <- (1-1/b)*diag(n)- (b+1/b-2)*a*P/(b-a)
  num <- ((1-a)^(pgamma*q/2))*(det(V))^(c/2)*(1-1/b)^(n*q/2)
  den <- (choose(p, pgamma))*(1 - a/b)^(pgamma*q/2)*(det(V+ t(Y)%*%zai%*%Y))^((n+c)/2)
  return(num/den)
}


s_53 <- function(X, y, temp = seq(1, 0.1, length.out = 20) , 
                 L = 20,J=20, burnin = burnin,g = g, M =20){
  
  #scaling
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(y)
  tau <- 1e-5
  
  
  # Initialize empty set of selected variables
  selected_vars <- 0
  model_a_store <- NULL
  model_r_store <- NULL
  accepted_models <- data.frame()
  
  # neighborhood function
  get_neighborhood <- function(selected_vars, X, y) {
    model <- lm(y ~ X[, selected_vars, drop=FALSE])  # Fit the model
    r <- residuals(model)
    h <- t(r) %*% X
    h_norms <- apply(h, 2, function(col) norm(col, type = '2'))
    scrs_vars <- order(h_norms, decreasing = TRUE)[1:M]
    setdiff(scrs_vars, selected_vars)  # Return only unselected variables
  }
  
  tic <- proc.time()[3]
  
  # choosing the initial model
  # Fitting a model and getting screened neighborhood
  model <- lm(y ~ 1) # k(1,1)
  r <- residuals(model)
  h <- t(r)%*%X
  h_norms <- apply(h, 2,function(col) norm(col,type = '2'))
  scrs_vars <- order(h_norms, decreasing = TRUE)[1:M]
  #model_a_store <- sapply(scrs_vars, function(s) c(p_addition = post(s), s = s))
  #selected_vars <- scrs_vars[which.max(model_a_store[1, ])]
  #accepted_models_list[[1]] <- list(post = max(model_a_store[1, ]), model = selected_vars, add = NA, rem = NA)
  
  
  remaining_vars <- setdiff(scrs_vars, selected_vars)
  
  # Only consider additive neighborhood for the initial model
  for (s in remaining_vars) {
    p_addition <- post(s)  # Compute the posterior probability
    model_a_store <- rbind(model_a_store, data.frame(p_addition,s))
  }
  selected_vars <- model_a_store$s[which.max(model_a_store$p_addition)]#map
  accepted_models <- data.frame(post = max(model_a_store$p_addition),
                                model = selected_vars,
                                add = NA,
                                rem = NA)
  iter <- 0
  iter.before.burn <- 0
  
  foo <- ceiling((L*J)*0.5)
  #starting the loop
  for (l in 1:L) {
    #print(paste0("Outer Loop Iteration l = ", l))
    for (i in 1:J){
      #print(paste0("  Inner Loop Iteration i = ", i))
      iter <- (l - 1) * J + i
      print(paste(l,i,iter, sep =";"))
      
      # Log current selected variables
      #cat("   Current Selected Vars: ", selected_vars, "\n")
      
      #for models with only one variable consider only additive nbd
      if(length(selected_vars) == 1){
        # Fitting a model and getting screened nbd
        remaining_vars <- get_neighborhood(selected_vars, X,Y)
        
        # Additive neighborhood and their posteriors
        model_a_store <- NULL
        for (s in remaining_vars){
          added_vars <- c(selected_vars, s )
          p_addition <-  post(added_vars)   # Compute the posterior probability
          model_a_store <- rbind(model_a_store, data.frame(p_addition, z = NA))
          model_a_store$z[nrow(model_a_store)] <- list(added_vars)
        }
        
        #model_a_store_list <- lapply(remaining_vars, function(s) {
        #  added_vars <- c(selected_vars, s)
        #  list(p_addition = post(added_vars), vars = added_vars)
        #})
        
        
        total_posta <- sum(model_a_store$p_addition)
        model_a_store$p_addition <- model_a_store$p_addition/total_posta
        accepted_models$add[nrow(accepted_models)] <- total_posta
        ad <- model_a_store[sample(nrow(model_a_store), 1, prob = (model_a_store$p_addition)^(1/temp[l])), ]#k+
        
        proposed_model <- ad #the additive model is proposed
        # Log proposed additive model
        #cat("    Proposed Additive Model (ad): ", unlist(ad$z), "\n")
        
        # accept-reject is not being used for this 
        #remove neighborhood posterior for this proposed model
        #if(length(unlist(ad$z))>1){
        #  total_postr_new <- 0
        #  for (s in unlist(ad$z)){
        #    reduced_vars <- setdiff(unlist(ad$z), s)
        #    total_postr_new <- total_postr_new + post(unlist(reduced_vars))
        #  }
        #}
        #logR <- log(accepted_models$add[nrow(accepted_models)]) - log(total_postr_new)#ratio for accept-reject
        #if(log(runif(1)) < logR){
        selected_vars <- unlist(proposed_model$z)
        accepted_models <- rbind(accepted_models, data.frame(post = proposed_model$p_addition,
                                                             model = NA,
                                                             add = NA,
                                                             rem = NA))
        accepted_models$model[nrow(accepted_models)] <- list(selected_vars)
        
        
        #}
      }else if(length(selected_vars) > 1){
        
        #for models having more than one variable
        
        # Fitting a model and getting screened nbd
        #model <- lm(y ~ X[, selected_vars]) 
        #r <- residuals(model)
        #h <- t(r)%*%X
        #h_norms <- apply(h, 2,function(col) norm(col,type = '2'))
        #scr_vars <- which(h_norms>mean(h_norms))
        remaining_vars <- get_neighborhood(selected_vars, X,Y)
        
        # Additive neighborhood and their posteriors
        model_a_store <- NULL
        for (s in remaining_vars){
          added_vars <- c(selected_vars, s )
          p_addition <-  post(added_vars)   # Compute the posterior probability
          model_a_store <- rbind(model_a_store, data.frame(p_addition, z = NA))
          model_a_store$z[nrow(model_a_store)] <- list(added_vars)
        }
        total_posta <- sum(model_a_store$p_addition)
        model_a_store$p_addition <- model_a_store$p_addition/total_posta
        accepted_models$add[nrow(accepted_models)] <- total_posta
        
        #to get posterior of any model use "post(unlist(model_a_store[1,2]))"
        #reduced neighborhood and their posteriors
        
        model_r_store <- NULL
        for (s in selected_vars){
          reduced_vars <- setdiff(selected_vars, s)
          p_reduced <- post(reduced_vars)   # Compute the posterior probability
          model_r_store <- rbind(model_r_store, data.frame(p_reduced, z = NA))
          model_r_store$z[nrow(model_r_store)] <- list(reduced_vars)
        }
        
        total_postr <- sum(model_r_store$p_reduced)
        model_r_store$p_reduced <- model_r_store$p_reduced/total_postr
        accepted_models$rem[nrow(accepted_models)] <- total_postr
        
        #sampling from additive and removed models
        ad <- model_a_store[sample(nrow(model_a_store), 1, 
                                   prob = (model_a_store$p_addition)^(1/temp[l])), ]#k+
        re <- model_r_store[sample(nrow(model_r_store), 1, 
                                   prob = (model_r_store$p_reduced)^(1/temp[l])), ]#k-
        
        
        # Log proposed models
        #cat("    Proposed Additive Model (ad): ", unlist(ad$z), "\n")
        #cat("    Proposed Removal Model (re): ", unlist(re$z), "\n")
        
        
        
        #sampling from additive and remove model
        if(sample(1:2,1, prob = c((ad[,1])^(1/temp[l]),(re[,1])^(1/temp[l]))) == 1){
          
          proposed_model <- ad #the additive model is proposed
          # accept-reject
          #remove neighborhood posterior for this proposed model
          total_postr_new <- 0
          for (s in unlist(proposed_model$z)){
            reduced_vars <- setdiff(unlist(proposed_model$z), s)
            total_postr_new <- total_postr_new + post(unlist(reduced_vars))
          }
          logR <- log(accepted_models$add[nrow(accepted_models)]) - log(total_postr_new)
          alpha <- min(1, exp(logR))
          #ratio for accept-reject
          if(runif(1) < alpha){
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
          # additive neighbourhood posterior for this proposed model
          model_a_store <- NULL
          remaining_vars <- get_neighborhood(unlist(proposed_model$z), X,Y)
          
          total_posta_new <-0
          for (s in remaining_vars){
            added_vars <- c(unlist(proposed_model$z), s )
            total_posta_new <- total_posta_new + post(added_vars)
          }
          logR <- log(accepted_models$rem[nrow(accepted_models)]) - log(total_posta_new)#ratio for accept-reject
          alpha <- min(1, exp(logR))
          if(runif(1) < alpha){
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
        # Log selected variables after update
        #cat("    Updated Selected Vars: ", selected_vars, "\n")
        
        # Stop if the number of selected variables reaches p-1
        if (length(selected_vars) == p - 1) {
          break
        }
      }
    }
    #print(paste0('iter.before.burn:',iter.before.burn))
    #print(accepted_models[,1:2])
    #print(proposed_model$z)#first loop ends
  }#outer loop ends
  
  toc <- proc.time()[3]
  a <- toc-tic
  print(a)
  #post-processing
  answer  <- accepted_models[-(1:iter.before.burn),1:2]
  post_burnin <- accepted_models[1:iter.before.burn,1:2]
  return(answer)
  
}

small <- list()
for (i in 1:1){
  answer_100 <- s_53(X,y=y, L=10,J=10,M=50)
  answer_200 <- s_53(X,y, L=10,J=20,M=20)
  answer_600 <- s_53(X,y, L=30,J=20,temp = seq(1, 0.1, length.out = 30))
  answer_1000<- s_53(X,y, L=10,J=100,temp = seq(1, 0.1, length.out = 100))
  
  small <- append(small,list(answer_100,answer_200,answer_600,answer_1000, constant) )
  names(small) <- c(100,200,600,1000, 'constant')
}
#for this we got gud result for g <- p^2/2 (= 5e+05)


##########################################################################
# Post MCMC


Bigg <- list(answer3_10,answer3_20,answer3_100)
names(Bigg) <- c(100,200,1000)
Bigg_40 <- Bigg

# the combination of variables that occurred the highest number of times
b <-  sort(table(sapply(answer$model, toString)), decreasing = T)
most_frequent_comb <- b[1]
most_frequent_comb <- unlist(names(most_frequent_comb))


#max posterior model
map_model <- answer[which.max(answer[,1]),]


### inclusion probability
prob = table(unlist(answer[,2]))/dim(answer)[1]
#median model
median_prob_model = as.numeric(names(prob[prob>0.5]))
