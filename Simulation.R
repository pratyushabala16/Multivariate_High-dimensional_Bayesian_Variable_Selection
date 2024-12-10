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
tau <-
g <- ##need to verify each time

# posterior probability
post<- function(selected_vars, c=3){
#to be updated
}


s_53 <- function(X, y, temp = seq(1, 0.1, length.out = 20) , 
                 L = 20,J=20, burnin = burnin,g = g, M =20){
#to be updated
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
