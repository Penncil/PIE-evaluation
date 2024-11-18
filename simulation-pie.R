## simulation under uniform prior distribution

#### SET UP parameters in the simulation #####
K <- 200  ##no.of data
n <- 15000
p.x <- 0.3
#p.x <- 0.05
alpha1 <- 0.65
alpha0 <- 0.99

## prevalence of outcome y given exposure/unexposure
#p.y.unex <- seq(1/100000, 5/100, length.out = 5)
p.y.unex <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
beta1 <- log(3)
beta0.values <- logit(p.y.unex)
beta0 <- beta0.values

delta <- c(0.05, 0.1, 0.15, 0.2)

p_len <- length(p.y.unex)
setting <- cbind(beta0[rep(1:p_len, each=4)], delta[rep(1:4,p_len)],delta[rep(1:4,p_len)])
setting <- rbind(setting, cbind(beta0[rep(1:p_len, each=2)],
                                matrix(c(0.05,0.15,0.15,0.05),2, byrow=TRUE)[rep(1:2,p_len),]))
colnames(setting) <- c('beta0', 'b','d')
setting <- data.frame(setting[order(setting[,1]),])

###### FUNC for GENERATING DATA ######

Data.gen <- function(n,p.x,beta0,beta1,alpha1,alpha0){
  x = rbinom(n, 1, p.x)   ## generate x from Bernoulli(0.5) distribution
  lin=beta0+beta1*x
  mu=exp(lin)/(1+exp(lin))
  y=rbinom(n,1,mu)
  sm=alpha1*(y==1)+(1-alpha0)*(y==0)
  s = rbinom(n,1,sm)
  list(s=s, x=x)
}


## run over 200 data generating from the same parameters
## Input: n, p.x, beta0, beta1, alpha1, alpha0, low.sen, up.sen, low.spe, up.spe
## Output: 200 rows of estimators & cov from inverse of the hessian

## task_id 1-200, for parallel computing on HPC
task_id <- 1

## fix the prior distribution of specificity, Uniform (0.95, 0.9999)
low2 <- 0.95; up2 <- 0.9999

for (ind in 1:dim(setting)[1]){
  parm <- setting[ind,]
  beta0 <- parm$beta0; b <- parm$b; d <- parm$d
  ## The prior distribution of sensitivity, Uniform (low1, up1)
  low1 <- alpha1 - b; up1 <- alpha1 + d;
  set.seed(211108+task_id*ind)
  
  data <- Data.gen(n, p.x, beta0, beta1, alpha1, alpha0)
  s <- data$s; x <- data$x
  
  ############ RUN THE MODEL ######################
  result <- tryCatch(NaiveVsPIE.uni.norm(s,x,adjust=1.4), 
                     warning = function(o) -1, error=function(o) -1)
  #################################################
  
  fail <- 0
  ## May fail due to numerical issues (see pie_func.R for detail)
  ## Regenerate data and refit the model if failed.
  ## If ALL the regenerated data failed, output -1
  while(is.numeric(result) & fail < 5){
    fail <- fail+1
    data <- Data.gen(n, p.x, beta0, beta1, alpha1, alpha0)
    s <- data$s; x <- data$x
    result <- tryCatch(NaiveVsPIE.uni.norm(s,x), warning = function(o) -1, error=function(o) -1)
  }
  if(fail > 4) output <- t(rep(-1,11))
  if(fail < 5) output <- t(c(result$out,fail))

  write.csv(output, paste0('/out-bias/s', ind, '/', task_id, ".csv"), row.names=FALSE)
}