# Load necessary packages/functions
library(BayesMallows)
library(Rcpp)
sourceCpp("Cpp/leapandshift_sourceCpp.cpp") # use the C++ leap-and-shift function

set.seed(17122021)

# Which dataset to simulate 
dataset="topK"
#dataset="rank_consistency"

# Generate data for simulation
N <- 10 # number of assessors
n <- 50 # total number of items
n_star <- 20 # number of items selected to be "relevant", i.e. will have the highest ranks
n_star_true <- n_star
alpha_true <- 2

# Tuning parmameters for the MCMC
M <- 1e4
leap_size = round(n_star/5)
L <- 1

# Simulate data 
if(dataset=="topK"){
  rho_true <- sample.int(n = n_star, size = n_star)
  data_red <- sample_mallows(rho_true, alpha_true, N, burnin = 1e4, thinning = 100) # ranks/data in A*
  A_star_true <- sort(sample.int(n = n, size = n_star)) # items in A*
  data <- matrix(NA, N, n) 
  data[,A_star_true] <- data_red
  true_rank <- A_star_true[sort(rho_true, index.return=TRUE)$ix]
  
  # Complete the data by randomly assigning to the other n-n* items the ranks from n*+1 to n
  for(j in 1:N)data[j,setdiff(1:n, A_star_true)] <- sample((n_star + 1):n, size = n-n_star)

}

if(dataset=="rank_consistency"){
  
  # Simulate data for the "relevant" items first
  rho_true <- sample.int(n = n_star, size = n_star)
  data_red <- sample_mallows(rho_true, alpha_true, N, burnin = 1e4, thinning = 100) # ranks/data in A*
  A_star_true <- sort(sample.int(n = n, size = n_star)) # items in A*
  true_rank <- A_star_true[sort(rho_true, index.return=TRUE)$ix]
  
  # The true rank order consistency
  rho_true_ordered <- sort(rho_true, index.return =TRUE)$ix
  
  # Randomly assign ranks while keeping order consistency
  for(j in 1:N)data_red[j,rho_true_ordered] <- sort(sample(1:n, size = n_star))
  
  # Complete the data by randomly assigning to the other n-n* items with ranks from 1,..,n
  data <- matrix(NA, N, n) 
  data[,A_star_true] <- data_red
  for(j in 1:N)data[j,setdiff(1:n, A_star_true)] <- sample(setdiff(1:n, data_red[j,]), size = n-n_star)
  
  
}

##############
### lowBMM ###
##############


# MCMC initialization: randomly generate rho0, alpha0 and A_star0
rho0 <- sample.int(n=n_star, size=n_star)
alpha0 <- alpha_true
A_star0 <- sort(sample.int(n=n, size=n_star))
prob_back <- prob_forw <- 0.5

# Initial data in dimension (N,n*)
init_data <- matrix(NA, N, n_star)
current_selection <- data[,A_star0] # select items from A*
for(i in 1:N)init_data[i, sort.int(current_selection[i,], index.return = TRUE)$ix] <- 1:n_star

# Proposal/old
prop_data <- old_data <- init_data
rho_prop <- rho_old <- rho0
A_star_prop <- A_star_old <- A_star0
alpha_prop <- alpha_old <- alpha0

# Compute the average rankings of all items to use for rho proposed
avg_ranks <- colMeans(data)

# MCMC
rho_mcmc <- matrix(NA, M, n)
rho_mcmc[1,A_star0] <- rho0
rho_prop_mcmc <- matrix(NA, M, n_star)
rho_prop_mcmc[1,] <- rho_prop
A_mcmc <- matrix(0, M, n)
A_mcmc[1,A_star0] <- 1
ACC_rho <- RATIO_rho <- ACC_A <- RATIO_A <- NULL

for(m in 2:M) {
  # MH step 1: update rho restricted on A*
  # 1a. Sample rank proposal through leap and shift (cpp function)
  rho_old <- rho_mcmc[m-1,A_star_old]
  tmp <- leap_and_shift(rho_proposal = rho_prop_mcmc[m-1,], indices = c(1:n_star), prob_backward = prob_back, prob_forward = prob_forw, rho = rho_old, leap_size = leap_size, reduce_indices = F)
  rho_prop <- tmp$rho_proposal[,]
  
  # 1b. Compute distances to current and proposed ranks
  dist_new <- abs(scale(old_data, rho_prop, scale = FALSE)) # footrule distance metric
  dist_old <- abs(scale(old_data, rho_old, scale = FALSE))
  rank_dist_sum <- sum(dist_new-dist_old)
  
  # 1c. Compute MH ratio and accept/reject
  prob_back <- tmp$prob_backward # probability backwards
  prob_forw <- tmp$prob_forward # probability forward
  prior_rho_old <- 1 # uniform prior
  prior_rho_prop <- 1
  C_rho <- (prob_back*prior_rho_old)/(prob_forw*prior_rho_prop) # leap and shift backwards/forwards probability factor (needs to be edited!!)
  ratio_rho <- min(1, C_rho*exp((-alpha_old/n_star)*rank_dist_sum))
  
  # Log version
  #ratio_rho <- log(prob_back)-log(prob_forw)-(alpha_old/n_star)*rank_dist_sum
  
  if(runif(1)<ratio_rho){
    rho_old <- rho_prop
    ACC_rho <- c(ACC_rho, 1)
  }else{
    ACC_rho <- c(ACC_rho, 0)
  }
  
  # MH step 2: update A*
  # 2a. Sample A* proposed by perturbing L items
  removed_items <- sample.int(n = n_star, size = L) # index of items to be removed in A*
  new_items <- sample(setdiff(1:n, A_star_old), size = L)
  A_star_prop <- sort(c(A_star_old[-removed_items], new_items))
  
  # 2b. Update data according to new set (ranks should go from in 1,..,n*)
  prop_selection <- data[,A_star_prop]
  for(i in 1:N)prop_data[i, sort.int(prop_selection[i,], index.return = TRUE)$ix] <- 1:n_star
  
  # 2c. Compute new corresponding rho_prop based on the items selected
  rho_prop_star <- rho_prop
  
  # Acceptance A*
  idx_match_old <- match(A_star_prop, A_star_old)[!is.na(match(A_star_prop, A_star_old))]
  idx_match_prop <- match(A_star_old, A_star_prop)[!is.na(match(A_star_old, A_star_prop))]
  rho_prop_star[idx_match_prop] <- rho_old[idx_match_old]
  idx <- match(new_items, A_star_prop) # index of the new item(s) in A*
  
  # Alternative: if we want to order the new items based on their avg data rankings
  #idx2 <- A_star_prop[idx]
  #idx3 <- sort(avg_ranks[idx2],index.return = TRUE)$ix # sorting new items based on ranking
  #rho_prop_star[idx[idx3]] <- rho_old[sort(removed_items)]
  
  # Alternative: randomly assign new items their ranking, not taking data into account
  rho_prop_star[idx] <- rho_old[removed_items]
  
  # 2d. Compute MH ratio and accept/reject
  ratio_A <- min(1, exp(-alpha_old/n_star*(sum(abs(scale(prop_data, rho_prop_star, scale = FALSE)))-
                                             sum(abs(scale(old_data, rho_old, scale = FALSE))))))
  if(runif(1)<ratio_A){
    A_star_old <- A_star_prop
    old_data <- prop_data
    rho_old <- rho_prop_star
    ACC_A <- c(ACC_A, 1)
  }else{
    ACC_A <- c(ACC_A, 0)
  }
  
  # Save current values
  rho_mcmc[m,A_star_old] <- rho_old
  rho_prop_mcmc[m,] <- rho_prop
  A_mcmc[m,A_star_old] <- 1
  RATIO_rho <- c(RATIO_rho, ratio_rho)
  RATIO_A <- c(RATIO_A, ratio_A)
  
  
  #if(m%%1000 == 0)print(paste('Iteration: ',m,sep=''))
}




#### RESULTS 

# Remove burnin 
burnin = M*0.05 # remove 5% 
rho_mcmc = rho_mcmc[c((burnin+1):M),]
M_burnin <- dim(rho_mcmc)[1] # M-burnin, how many samples we have for post-processing

# Marginal posterior mean of rho 
avg_rho_ranks = colSums(rho_mcmc, na.rm = TRUE)/colSums(!is.na(rho_mcmc))
rankedItems <- sort(avg_rho_ranks, index.return=T, na.last = T)$ix

# Proportion of correct items in the set
tmp <- match(rankedItems[1:n_star], true_rank)
no_correct_items <- length(tmp[!is.na(tmp)])
no_misclassified <- n_star - no_correct_items
prop_correct_items <- no_correct_items/n_star 

# Save results to res: footrule distance, proportion of correct items selected
res <- data.frame(index=numeric(1))
res$distance[[1]] <- abs(true_rank-rankedItems[1:n_star])
res$prop_correct <- prop_correct_items

# Save computational times
res$comp_time <- as.character(comp_time_lowbmm)

# Compute the additional two measures from Zhu et al. (2021):
# Recovery distance: R_rho = d(rho, rho_true) + no_misclasified*0.5*(n + n* + 1)
res$recovery_dist <- rank_distance(true_rank, rankedItems[1:n_star], metric="kendall") + no_misclassified*0.5*(n + n_star + 1)
# Coverage distance: c_rho = (n* - no_misclassified)/n*
res$coverage <- (n_star-no_misclassified)/n_star # same as prop_correct

# Save method name
res$model <- "lowBMM"

# Save results for postprocessing and plotting
#save(res, n_star, n_star_true, L, leap_size, alpha_true, M, n, N, file=sprintf('compare_methodsBIG_%s_alphafixed%g_nstartrue%g_nstar%g_L%d_leap%d_n%d_N%d_M%d_reps%d.RData', dataset, alpha_true, n_star_true, n_star, L, leap_size, n, N, M, R))


