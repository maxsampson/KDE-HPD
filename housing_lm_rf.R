rm(list = ls())

library(hdrcde)
realestate <- read.delim("realestate.txt") #https://online.stat.psu.edu/stat501/lesson/8/8.9
dim(realestate)
mod <- lm(SalePrice ~ SqFeet + Air + SqFeet:Air, data = realestate)
pdf("resid_home.pdf", width = 8, height = 4)
plot(mod$residuals ~ realestate$SqFeet, ylab = "Residuals", xlab = "Square Feet (in thousands)")
dev.off()
score_fun <- function(y_hat, y_true, sig){ #for upper one-sided, add for upper, subtract for lower
    (y_true - y_hat) / sig
}
quant_fun <- function(dat, alpha, n_size, direction){
    o <- order(dat)
    if(all.equal(alpha, 1) == TRUE) {quantile(dat, 1)}
    else if(all.equal(alpha, 0) == TRUE) {quantile(dat, 0)}
    else if(direction == "low"){dat[o][ceiling((alpha) * (n_size + 1))]}
    else if(direction != "low"){dat[o][ceiling((alpha) * (n_size + 1)) - 1]}
}

## KDE-HPD
perms <- 200
cov <- matrix(NA, nrow = perms, ncol = 221) 
len <- matrix(NA, nrow = perms, ncol = 221) 
cov_m <- matrix(NA, nrow = perms, ncol = 221)
cov_f <- matrix(NA, nrow = perms, ncol = 221)
cov_large <- matrix(NA, nrow = perms, ncol = 221)

cov_parm <- matrix(NA, nrow = perms, ncol = 221)
len_parm <- matrix(NA, nrow = perms, ncol = 221) 
cov_parm_m <- c()
cov_parm_f <- c()

set.seed(1)
for(i in seq(perms)){
    indicies <- sample(521, 521, FALSE)
    dat <- realestate[indicies, ]
    dat_train <- dat[1:60, ] #500 for mean model
    dat_hetero <- dat[61:200, ] #500 for heteroskedastic 
    dat_cal <- dat[201:300, ] #100 calibration
    dat_out <- dat[301:521, ] #out of sample test
    
    mod_train <- lm(SalePrice ~ SqFeet + Air + SqFeet:Air, data = dat_train)
    
    hetero_test <- predict(mod_train, newdata = dat_hetero)
    sigma_test <- abs(hetero_test - dat_hetero$SalePrice)
    dat_hetero <- cbind(dat_hetero, sigma_test)
    
    mod_sig <- ranger::ranger(sigma_test ~ SqFeet + Air, data = dat_hetero)
    sig_cal <- predict(mod_sig, data = dat_cal)$predictions
    sig_out <- predict(mod_sig, data = dat_out)$predictions

    pred_cal <- predict(mod_train, newdata = dat_cal)
    
    upper_scores <- score_fun(pred_cal, dat_cal$SalePrice, (sig_cal))
    lower_scores <- -1 * upper_scores
    
    
    dens <- density(lower_scores, kernel = "gaussian", from = range(lower_scores)[1], to = range(lower_scores)[2], n = 500, adjust = 100 ^ (-1/3) / 100 ^ (-1/5))
    hpd_result <- hdr(lower_scores, prob = 90, den = dens)
    
    cdf <- ecdf(lower_scores)
    if(length(hpd_result$hdr) > 4){
        
        lower_alpha1 <- cdf(max(hpd_result$hdr)) - cdf(hpd_result$hdr[3]) + cdf(hpd_result$hdr[2]) - 0.9
        lower_alpha2 <- cdf(hpd_result$hdr[3])
        
        upper_alpha1 <- cdf(hpd_result$hdr[2]) 
        upper_alpha2 <- cdf(max(hpd_result$hdr)) 
        
        upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, length(lower_scores), "upper")
        lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, length(lower_scores), "low")
        
        upper_error2 <- quant_fun(lower_scores, 1-upper_alpha2, length(lower_scores), "upper")
        lower_error2 <- quant_fun(lower_scores, 1-lower_alpha2, length(lower_scores), "low")
        
        pred_out <- predict(mod_train, newdata = dat_out)
        
        upper_out1 <- pred_out - upper_error1 * abs(sig_out)
        lower_out1 <- pred_out - lower_error1 * abs(sig_out)
        upper_out2 <- pred_out - upper_error2 * abs(sig_out)
        lower_out2 <- pred_out - lower_error2 * abs(sig_out)
        cov[i, ] <- (dat_out$SalePrice > lower_out1 & dat_out$SalePrice < upper_out1 | dat_out$SalePrice > lower_out2 & dat_out$SalePrice < upper_out2) 
        len[i, ] <- (upper_out1 - lower_out1 + upper_out2 - lower_out2)
    } else{
        lower_alpha1 <- cdf(max(hpd_result$hdr)) - 0.9
        
        upper_alpha1 <- cdf(max(hpd_result$hdr)) 
        
        upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, length(lower_scores), "upper")
        lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, length(lower_scores), "low")
        
        pred_out <- predict(mod_train, newdata = dat_out)
        
        upper_out1 <- (pred_out - upper_error1 * abs(sig_out))
        lower_out1 <- (pred_out - lower_error1 * abs(sig_out))
        temp <- mean(upper_out1 - lower_out1)
        cov[i, ] <- (dat_out$SalePrice > lower_out1 & dat_out$SalePrice < upper_out1) 
        len[i, ] <- upper_out1 - lower_out1
    }
    cov_m[i, ] <- c(cov[i, dat_out$Air == 1], rep(NA, 221 - sum(dat_out$Air == 1)))
    cov_f[i, ] <- c(cov[i, dat_out$Air == 0], rep(NA, 221 - sum(dat_out$Air == 0)))
    cov_large[i, ] <- c(cov[i, dat_out$SalePrice > 335], rep(NA, 221 - sum(dat_out$SalePrice > 335)))
    
    param_dat <- rbind(dat_train, dat_hetero[, 1:12], dat_cal)
    linear_mod <- lm(SalePrice ~ SqFeet + Air + SqFeet:Air, data = param_dat)
    
    param_interval <- predict(linear_mod, newdata = dat_out, interval = "predict", level = 0.90)
    cov_parm[i, ] <- dat_out$SalePrice > param_interval[, 2] & dat_out$SalePrice < param_interval[, 3] 
    len_parm[i, ] <- mean(param_interval[, 3] - param_interval[, 2])
    


    
    cov_parm_m <- c(cov_parm_m, cov_parm[i, dat_out$Air == 1])
    cov_parm_f <- c(cov_parm_f, cov_parm[i, dat_out$Air == 0])
}
mean(cov)
mean(len)
mean(cov_m, na.rm = TRUE)
mean(cov_f, na.rm = TRUE)
mean(cov_large, na.rm = TRUE)
median(apply(len, 1, median))


sd(rowMeans(cov)) / sqrt(200)
sd(rowMeans(len)) / sqrt(200)
sd(apply(len, 1, median)) / sqrt(200)
sd(rowMeans(cov_m, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
sd(rowMeans(cov_f, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
sd(rowMeans(cov_large, na.rm = TRUE), na.rm = TRUE) / sqrt(200)


## HPD-Split

cov2 <- matrix(NA, nrow = perms, ncol = 221) 
len2 <- matrix(NA, nrow = perms, ncol = 221) 
cov2_m <- matrix(NA, nrow = perms, ncol = 221)
cov2_f <- matrix(NA, nrow = perms, ncol = 221)
cov2_large <- matrix(NA, nrow = perms, ncol = 221)
library(predictionBands)
library(FlexCoDE)
set.seed(1)
for(i in seq(perms)){
    indicies <- sample(521, 521, FALSE)
    dat <- realestate[indicies, ]
    dat_train <- dat[1:300, ] #500 for mean model
    dat_out <- dat[301:521, ] #out of sample test

    fit <- fit_predictionBands(cbind(as.matrix(dat_train$SqFeet), as.matrix(dat_train$Air)), dat_train$SalePrice)
    
    bands <- predict(fit, cbind(as.matrix(dat_out$SqFeet), as.matrix(dat_out$Air)), type = "hpd", alpha = 0.1)
    
    temp_cov <- rep(NA, length(dat_out$SalePrice))
    temp_len <- rep(NA, length(dat_out$SalePrice))
    for(j in seq(length(dat_out$SalePrice))){
        keep <- c(FALSE, unlist(bands$prediction_bands_which_belong[j]), FALSE)
        grids <- c(min(dat$SalePrice), bands$y_grid, max(dat$SalePrice))
        interval <- grids[keep]
        if(sum(abs(diff(keep))) > 2 & sum(abs(diff(keep))) < 5){
            index <- which(diff(keep) != 0) ##assuming that the first change is False to True
            temp_cov[j] <- dat_out$SalePrice[j] > grids[index[1]] & dat_out$SalePrice[j] < grids[index[2]] | 
                dat_out$SalePrice[j] > grids[index[3]] & dat_out$SalePrice[j] < grids[index[4]]
            temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] 
        } else if(sum(abs(diff(keep))) > 5){
            index <- which(diff(keep) != 0) ##assuming that the first change is False to True
            temp_cov[j] <- dat_out$SalePrice[j] > grids[index[1]] & dat_out$SalePrice[j] < grids[index[2]] | 
                dat_out$SalePrice[j] > grids[index[3]] & dat_out$SalePrice[j] < grids[index[4]] |
                dat_out$SalePrice[j] > grids[index[5]] & dat_out$SalePrice[j] < grids[max(index)]
            temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] + grids[max(index)] - grids[index[5]]
        } else{
            temp_cov[j] <- dat_out$SalePrice[j] > min(interval) & dat_out$SalePrice[j] < max(interval)
            temp_len[j] <- max(interval) - min(interval)
        }
    }
    
    cov2[i, ] <- temp_cov 
    len2[i, ] <- (temp_len)
    
    cov2_m[i, ] <- c(cov2[i, dat_out$Air == 1], rep(NA, 221 - sum(dat_out$Air == 1)))
    cov2_f[i, ] <- c(cov2[i, dat_out$Air == 0], rep(NA, 221 - sum(dat_out$Air == 0)))
    cov2_large[i, ] <- c(cov2[i, dat_out$SalePrice > 335], rep(NA, 221 - sum(dat_out$SalePrice > 335)))
    
    print(i)
}
mean(cov2)
mean(len2)
median(apply(len2, 1, median))
mean(cov2_m, na.rm = TRUE)
mean(cov2_f, na.rm = TRUE)
mean(cov2_large, na.rm = TRUE)

sd(rowMeans(cov2)) / sqrt(200)
sd(rowMeans(len2)) / sqrt(200)
sd(apply(len2, 1, median)) / sqrt(200)
sd(rowMeans(cov2_m, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
sd(rowMeans(cov2_f, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
sd(rowMeans(cov2_large, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
median(len2)
## CQR
alpha.sig <- 0.1
cqr <- function(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig){
    
    XXX <- rbind(X1,X.test)
    YYY <- c(Y1,Y.test)
    
    quant_mod <- quantregRanger(y ~ . , data = data.frame(y = Y0, x = X0))
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

library(quantreg)
library(quantregRanger)

cov_mat_cqr <- matrix(NA, nrow = perms, ncol = 221) 
len_mat_cqr <- matrix(NA, nrow = perms, ncol = 221) 
cov_cqr_m <- matrix(NA, nrow = perms, ncol = 221)
cov_cqr_f <- matrix(NA, nrow = perms, ncol = 221)
cov_cqr_large <- matrix(NA, nrow = perms, ncol = 221)

set.seed(1)
for (i in seq(perms)){
    indicies <- sample(521, 521, FALSE)
    dat <- realestate[indicies, ]
    dat_train <- dat[1:200, ] #500 for mean model
    dat_cal <- dat[201:300, ] #1000 calibration
    dat_out <- dat[301:521, ] #out of sample test
    
    # Define training and holdout samples
    
    Y.test  <- dat_out$SalePrice
    X.test  <- cbind(dat_out$Air, dat_out$SqFeet)
    X0 <- cbind(dat_train$Air, dat_train$SqFeet) ## train
    Y0 <- dat_train$SalePrice ## train
    X1 <- cbind(dat_cal$Air, dat_cal$SqFeet) ## calibrate
    Y1 <- dat_cal$SalePrice ## calibrate
    
    # cqr
    res.opt   <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
    
    cov_mat_cqr[i, ]   <- res.opt$cov.opt
    len_mat_cqr[i, ]  <- res.opt$leng.opt
    
    cov_cqr_m[i, ] <- c(cov[i, dat_out$Air == 1], rep(NA, 221 - sum(dat_out$Air == 1)))
    cov_cqr_f[i, ] <- c(cov[i, dat_out$Air == 0], rep(NA, 221 - sum(dat_out$Air == 0)))
    cov_cqr_large[i, ] <- c(cov[i, dat_out$SalePrice > 335], rep(NA, 221 - sum(dat_out$SalePrice > 335)))
    
    print(i)
}

# Average length
mean(cov_mat_cqr)
mean(len_mat_cqr)
median(apply(len_mat_cqr, 1, median))
mean(cov_cqr_m, na.rm = TRUE)
mean(cov_cqr_f, na.rm = TRUE)
mean(cov_cqr_large, na.rm = TRUE)

sd(rowMeans(cov_mat_cqr)) / sqrt(200)
sd(rowMeans(len_mat_cqr)) / sqrt(200)
sd(apply(len_mat_cqr, 1, median)) / sqrt(200)
sd(rowMeans(cov_cqr_m, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
sd(rowMeans(cov_cqr_f, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
sd(rowMeans(cov_cqr_large, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
median(len_mat_cqr)


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
    lower <- NULL; upper <- NULL
    
    for (t in ind.test){
        ci.grid <- abs(taus - bhat[t]-(1-alpha.sig)/2)
        ci <- Q.yx[t,(ci.grid <= threshold)]
        ub <- max(ci)
        lb <- min(ci)
        leng.opt <- c(leng.opt,ub-lb)
        lower <- c(lower, lb)
        upper <- c(upper, ub)
        
    }
    
    leng.opt[which(leng.opt==-Inf)] <- NA
    
    return(list(cov.opt=cov.opt,leng.opt=leng.opt,lower=lower,upper=upper))  
    
}

cov_mat_dcp <- matrix(NA, nrow = perms, ncol = 221) 
len_mat_dcp <- matrix(NA, nrow = perms, ncol = 221)
cov_dcp_m <- matrix(NA, nrow = perms, ncol = 221)
cov_dcp_f <- matrix(NA, nrow = perms, ncol = 221)
cov_large_dcp <- matrix(NA, nrow = perms, ncol = 221)

set.seed(1)
for (i in seq(perms)){
    indicies <- sample(521, 521, FALSE)
    dat <- realestate[indicies, ]
    dat_train <- dat[1:200, ] #500 for mean model
    dat_cal <- dat[201:300, ] #1000 calibration
    dat_out <- dat[301:521, ] #out of sample test
    
    # Define training and holdout samples
    
    Y.test  <- dat_out$SalePrice
    X.test  <- cbind(dat_out$Air, dat_out$SqFeet)
    X0 <- cbind(dat_train$Air, dat_train$SqFeet) ## train
    Y0 <- dat_train$SalePrice ## train
    X1 <- cbind(dat_cal$Air, dat_cal$SqFeet) ## calibrate
    Y1 <- dat_cal$SalePrice ## calibrate
    
    # cqr
    taus    <- seq(0.001,0.999,length=200)
    
    res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
    
    cov_mat_dcp[i, ]   <- res.opt$cov.opt
    len_mat_dcp[i, ]  <- res.opt$leng.opt
    
    cov_dcp_m[i, ] <- c(cov_mat_dcp[i, dat_out$Air == 1], rep(NA, 221 - sum(dat_out$Air == 1))) ##so that we can compute SD
    cov_dcp_f[i, ] <- c(cov_mat_dcp[i, dat_out$Air == 0], rep(NA, 221 - sum(dat_out$Air == 0))) ##so that we can compute SD
    cov_large_dcp[i, ] <- c(cov_mat_dcp[i, dat_out$SalePrice > 335], rep(NA, 221 - sum(dat_out$SalePrice > 335))) ##so that we can compute SD
    
    print(i)
}
mean(cov_mat_dcp)
mean(len_mat_dcp)
median(apply(len_mat_dcp, 1, median))
mean(cov_dcp_m, na.rm = TRUE)
mean(cov_dcp_f, na.rm = TRUE)
mean(cov_large_dcp, na.rm = TRUE)

sd(rowMeans(cov_mat_dcp)) / sqrt(200)
sd(rowMeans(len_mat_dcp)) / sqrt(200)
sd(apply(len_mat_dcp, 1, median)) / sqrt(200)
sd(rowMeans(cov_dcp_m, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
sd(rowMeans(cov_dcp_f, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
sd(rowMeans(cov_large_dcp, na.rm = TRUE), na.rm = TRUE) / sqrt(200)
median(len_mat_dcp)


#save.image(file = "housing_lm_rf.RDATA")


