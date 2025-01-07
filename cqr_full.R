###################################################################################
# Paper: Distributional conformal prediction
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# Empirical illustration: predicting daily returns
# DISCLAIMER: This software is provided "as is" without warranty of 
# any kind, expressed or implied.
###################################################################################

# Preliminaries
rm(list = ls())
cqr <- function(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig){
    
    XXX <- rbind(X1,X.test)
    YYY <- c(Y1,Y.test)
    
    quant_mod <- quantregRanger(y ~ x, data = data.frame(y = Y0, x = X0))
    yhat <- predict(quant_mod, data = data.frame(x = XXX), quantiles = c(alpha.sig / 2, 1 - alpha.sig / 2))
    score_vec <- pmax(yhat[1:length(Y1), 1] - Y1, Y1 - yhat[1:length(Y1), 2])
    
    k <- ceiling((1 - alpha.sig) * (1 + length(Y1)))
    final_adj <- sort(score_vec)[k]
    
    lb.o <- yhat[(length(Y1) + 1):length(YYY), 1] - final_adj
    ub.o <- yhat[(length(Y1) + 1):length(YYY), 2] + final_adj
    

    cov.o <- (Y.test <= ub.o & Y.test >= lb.o)
    
    leng.o <- ub.o - lb.o
    
    
    return(list(cov.opt=cov.o,leng.opt=leng.o))  
    
}

# Working directory

# Packages
library(quantregRanger)

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
        X.test  <- cbind(x[1001:1050])
        X0 <- cbind(x[1:500]) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000]) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # cqr
        res.opt   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat1   <- rbind(cov.mat1,cov.mat.temp)
        leng.mat1  <- rbind(leng.mat1,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat1)
colMeans(leng.mat1)

# unimodal and symmetric

cov.mat2 <- leng.mat2 <- NULL
set.seed(6)
time_2_dcp <- system.time(
    for (i in seq(N)){
        x <- runif(1050, -5, 5)
        error <- rgamma(1050, shape = 7.5, rate = 1 )
        y <- 5 + 2 * x + error
        
        # Define training and holdout samples
        
        Y.test  <- y[1001:1050]
        X.test  <- cbind(x[1001:1050])
        X0 <- cbind(x[1:500]) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000]) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # cqr
        res.opt   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat2  <- rbind(cov.mat2,cov.mat.temp)
        leng.mat2  <- rbind(leng.mat2,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat2)
colMeans(leng.mat2)

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
        X.test  <- cbind(x[1001:1050])
        X0 <- cbind(x[1:500]) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000]) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # cqr

        res.opt   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat3  <- rbind(cov.mat3,cov.mat.temp)
        leng.mat3  <- rbind(leng.mat3,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat3)
colMeans(leng.mat3)

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
        X.test  <- cbind(x[1001:1050])
        X0 <- cbind(x[1:500]) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000]) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # cqr
        res.opt   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat4  <- rbind(cov.mat4,cov.mat.temp)
        leng.mat4  <- rbind(leng.mat4,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat4)
colMeans(leng.mat4)

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
        X.test  <- cbind(x[1001:1050])
        X0 <- cbind(x[1:500]) ## train
        Y0 <- y[1:500] ## train
        X1 <- cbind(x[501:1000]) ## calibrate
        Y1 <- y[501:1000] ## calibrate
        
        # cqr
        res.opt   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
        
        cov.mat.temp  <- mean(res.opt$cov.opt)
        leng.mat.temp <- mean(res.opt$leng.opt)
        cov.mat5  <- rbind(cov.mat5,cov.mat.temp)
        leng.mat5  <- rbind(leng.mat5,leng.mat.temp)
        print(i)
    })

# Average length
colMeans(cov.mat5)
colMeans(leng.mat5)

#save.image(file = "cqr_simulations_run.RData")

