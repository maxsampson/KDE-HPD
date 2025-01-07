###################################################################################
# Paper: Distributional conformal prediction
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# Empirical illustration: predicting daily returns
# DISCLAIMER: This software is provided "as is" without warranty of 
# any kind, expressed or implied.
###################################################################################

# Preliminaries
rm(list = ls())
dcp.opt <- function(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig){
    
    XXX <- rbind(X1,X.test)
    YYY <- c(Y1,Y.test)
    
    beta.qr <- matrix(NA,dim(X0)[2]+1,length(taus))
    for (t in 1:length(taus)){
        beta.qr[,t] <- rq.fit.br(cbind(1,X0),Y0,tau=taus[t])$coefficients ## coefficients
    }
    tQ.yx   <- cbind(1,XXX)%*%beta.qr ## yhat, each column is a new tau
    Q.yx    <- t(apply(tQ.yx,1,FUN=sort))
    u.hat <- rowMeans((Q.yx <= matrix(YYY,length(YYY),dim(tQ.yx)[2])))
    
    bhat <- rep(NA,length(YYY))
    b.grid <- taus[taus<=alpha.sig]
    for (t in 1:length(YYY)){
        leng <- rep(NA,length(b.grid))
        leng.test <- rep(NA,length(b.grid))
        for (b in 1:length(b.grid)){
            Q.yx.u <- approx(x=taus,y=Q.yx[t,],xout=(b.grid[b]+1-alpha.sig),rule=2)$y
            leng[b] <- Q.yx.u -Q.yx[t,b]
        }
        bhat[t] <- b.grid[which.min(leng)]
    }
    
    ind.test <- (length(Y1)+1):length(YYY)
    
    cs.opt <- abs(u.hat-bhat-(1-alpha.sig)/2)
    
    k           <- ceiling((1-alpha.sig)*(1+length(Y1)))
    threshold   <- sort(cs.opt[-ind.test])[k]
    
    cov.opt   <- (cs.opt[ind.test] <= threshold)
    
    leng.opt <- NULL
    for (t in ind.test){
        ci.grid <- abs(taus - bhat[t]-(1-alpha.sig)/2)
        ci <- Q.yx[t,(ci.grid <= threshold)]
        ub <- max(ci)
        lb <- min(ci)
        leng.opt <- c(leng.opt,ub-lb)
    }
    
    leng.opt[which(leng.opt==-Inf)] <- NA
    
    return(list(cov.opt=cov.opt,leng.opt=leng.opt))  
    
}

# Working directory

# Packages
library(quantreg)

# Functions
#source("functions_final.R")

###################################################################################
# Analysis
###################################################################################

# Miscoverage level
alpha.sig <- 0.1 
N <- 1e3

# unimodal and symmetric

cov.mat1 <- leng.mat1 <- NULL
set.seed(6)
time_1_dcp <- system.time(
for (i in seq(N)){
    x <- runif(1050, -5, 5)
    error <- rnorm(1050, 0, 1)
    y <- 5 + 2 * x + error
    
    # Define training and holdout samples
    
    Y.test  <- y[1001:1050]
    X.test  <- cbind(x[1001:1050], x[1001:1050]^2)
    X0 <- cbind(x[1:500], x[1:500]^2) ## train
    Y0 <- y[1:500] ## train
    X1 <- cbind(x[501:1000], x[501:1000]^2) ## calibrate
    Y1 <- y[501:1000] ## calibrate
    
    # DCP opt
    taus    <- seq(0.001,0.999,length=200)

    res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)

    cov.mat.temp  <- mean(res.opt$cov.opt)
    leng.mat.temp <- mean(res.opt$leng.opt)
    cov.mat1   <- rbind(cov.mat1,cov.mat.temp)
    leng.mat1  <- rbind(leng.mat1,leng.mat.temp)
    print(i)
})

# Average length
colMeans(cov.mat1)
colMeans(leng.mat1)
sd(cov.mat1) / sqrt(N)
sd(leng.mat1) / sqrt(N)

# unimodal and skewed

cov.mat2 <- leng.mat2 <- NULL
set.seed(6)
time_2_dcp <- system.time(
    for (i in seq(N)){
        x <- runif(1050, -5, 5)
        error <- rgamma(1050, shape = 7.5, rate = 1 )
        y <- 5 + 2 * x + error
        
        # Define training and holdout samples
        
        Y.test  <- y[1001:1050]
        X.test  <- cbind(x[1001:1050], x[1001:1050]^2)
        X0 <- cbind(x[1:500], x[1:500]^2) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000], x[501:1000]^2) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # DCP opt
        taus    <- seq(0.001,0.999,length=200)

        res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat2  <- rbind(cov.mat2,cov.mat.temp)
        leng.mat2  <- rbind(leng.mat2,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat2)
colMeans(leng.mat2)
sd(cov.mat2) / sqrt(N)
sd(leng.mat2) / sqrt(N)

# bimodal

cov.mat3 <- leng.mat3 <- NULL
set.seed(6)
time_3_dcp <- system.time(
    for (i in seq(N)){
        x <- runif(1050, -5, 5)
        index <- sample(1050, 1050, replace = FALSE)
        error1 <- rnorm(1050/2, -6)
        error2 <- rnorm(1050/2, 6)
        error <- c(error1, error2)[index]
        y <- 5 + 2 * x + error
        
        # Define training and holdout samples
        
        Y.test  <- y[1001:1050]
        X.test  <- cbind(x[1001:1050], x[1001:1050]^2)
        X0 <- cbind(x[1:500], x[1:500]^2) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000], x[501:1000]^2) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # DCP opt
        taus    <- seq(0.001,0.999,length=200)

        res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat3  <- rbind(cov.mat3,cov.mat.temp)
        leng.mat3  <- rbind(leng.mat3,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat3)
colMeans(leng.mat3)
sd(cov.mat3) / sqrt(N)
sd(leng.mat3) / sqrt(N)
# heteroskedastic

cov.mat4 <- leng.mat4 <- NULL
set.seed(6)
time_4_dcp <- system.time(
    for (i in seq(N)){
        x <- runif(1050, -5, 5)
        error <- rgamma(1050, shape = 1 + 2 * abs(x), rate = 1 + 2 * abs(x))
        y <- 5 + 2 * x + error
        
        # Define training and holdout samples
        
        Y.test  <- y[1001:1050]
        X.test  <- cbind(x[1001:1050], x[1001:1050]^2)
        X0 <- cbind(x[1:500], x[1:500]^2) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000], x[501:1000]^2) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # DCP opt
        taus    <- seq(0.001,0.999,length=200)

        res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat4  <- rbind(cov.mat4,cov.mat.temp)
        leng.mat4  <- rbind(leng.mat4,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat4)
colMeans(leng.mat4)
sd(cov.mat4) / sqrt(N)
sd(leng.mat4) / sqrt(N)
# bowtie

cov.mat5 <- leng.mat5 <- NULL
set.seed(6)
time_5_dcp <- system.time(
    for (i in seq(N)){
        x <- runif(1050, -5, 5)
        error <- error <- rnorm(1050, 0, abs(x))
        y <- 5 + 2 * x + error
        
        # Define training and holdout samples
        
        Y.test  <- y[1001:1050]
        X.test  <- cbind(x[1001:1050], x[1001:1050]^2)
        X0 <- cbind(x[1:500], x[1:500]^2) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000], x[501:1000]^2) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # DCP opt
        taus    <- seq(0.001,0.999,length=200)

        res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat5  <- rbind(cov.mat5,cov.mat.temp)
        leng.mat5  <- rbind(leng.mat5,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat5)
colMeans(leng.mat5)
sd(cov.mat5) / sqrt(N)
sd(leng.mat5) / sqrt(N)
save.image(file = "dcp_simulations_run.RData")

