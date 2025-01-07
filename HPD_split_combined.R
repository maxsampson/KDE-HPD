rm(list = ls())

library(predictionBands)
library(FlexCoDE)

N <- 1e3

## here we look at a standard normal error
cov2_uni_norm <- matrix(NA, nrow = N, ncol = 50) 
len2_uni_norm <- matrix(NA, nrow = N, ncol = 50) 

set.seed(6)
time_theirs_norm <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500
        x1 <- runif(n, -5, 5)
        error <- rnorm(n)
        y <- 5 + 2 * x1 + error
        
        dat_train <- data.frame(y = y[1:1000], x1 = x1[1:1000])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x1 = x1[(ntest + ntrain + 1):n])
        
        fit <- fit_predictionBands(dat_train$x1, dat_train$y)
        
        bands <- predict(fit, as.matrix(dat_out$x1), type = "hpd", alpha = 0.1)
        
        temp_cov <- rep(NA, length(dat_out$y))
        temp_len <- rep(NA, length(dat_out$y))
        for(j in seq(length(dat_out$y))){
            keep <- c(FALSE, unlist(bands$prediction_bands_which_belong[j]), FALSE)
            grids <- c(min(y), bands$y_grid, max(y))
            interval <- grids[keep]
            if(sum(abs(diff(keep))) > 2 & sum(abs(diff(keep))) < 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] 
            } else if(sum(abs(diff(keep))) > 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]] |
                    dat_out$y[j] > grids[index[5]] & dat_out$y[j] < grids[max(index)]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] + grids[max(index)] - grids[index[5]]
            } else{
                temp_cov[j] <- dat_out$y[j] > min(interval) & dat_out$y[j] < max(interval)
                temp_len[j] <- max(interval) - min(interval)
            }
        }
        
        cov2_uni_norm[i, ] <- (temp_cov) #quantile regression and upper/lower adjustments
        len2_uni_norm[i, ] <- (temp_len)
        
        print(i)
    }
)

mean(cov2_uni_norm)
sd(rowMeans(cov2_uni_norm)) / sqrt(N)
mean(len2_uni_norm)
sd(rowMeans(len2_uni_norm)) / sqrt(N)

## here we look at a skewed error independent of X (unimodal and skewed)

cov2_uni_ind <- matrix(NA, nrow = N, ncol = 50) 
len2_uni_ind <- matrix(NA, nrow = N, ncol = 50) 

set.seed(6)
time_theirs_ind <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500
        x1 <- runif(n, -5, 5)
        error <- rgamma(n, shape = 7.5, rate = 1 )
        y <- 5 + 2 * x1 + error
        
        dat_train <- data.frame(y = y[1:1000], x1 = x1[1:1000])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x1 = x1[(ntest + ntrain + 1):n])
        
        fit <- fit_predictionBands(dat_train$x1, dat_train$y)
        
        bands <- predict(fit, as.matrix(dat_out$x1), type = "hpd", alpha = 0.1)
        
        temp_cov <- rep(NA, length(dat_out$y))
        temp_len <- rep(NA, length(dat_out$y))
        for(j in seq(length(dat_out$y))){
            keep <- c(FALSE, unlist(bands$prediction_bands_which_belong[j]), FALSE)
            grids <- c(min(y), bands$y_grid, max(y))
            interval <- grids[keep]
            if(sum(abs(diff(keep))) > 2 & sum(abs(diff(keep))) < 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] 
            } else if(sum(abs(diff(keep))) > 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]] |
                    dat_out$y[j] > grids[index[5]] & dat_out$y[j] < grids[max(index)]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] + grids[max(index)] - grids[index[5]]
            } else{
                temp_cov[j] <- dat_out$y[j] > min(interval) & dat_out$y[j] < max(interval)
                temp_len[j] <- max(interval) - min(interval)
            }
        }
        
        cov2_uni_ind[i, ] <- (temp_cov) #quantile regression and upper/lower adjustments
        len2_uni_ind[i, ] <- (temp_len)
        
        print(i)
    }
)

mean(cov2_uni_ind)
sd(rowMeans(cov2_uni_ind)) / sqrt(N)
mean(len2_uni_ind)
sd(rowMeans(len2_uni_ind)) / sqrt(N)

## here we look at a bimodal error term, independent of X (bimodal)

cov2 <- matrix(NA, nrow = N, ncol = 50) 
len2 <- matrix(NA, nrow = N, ncol = 50) 

set.seed(6)
time_theirs <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500
        x1 <- runif(n, -5, 5)
        error1 <- rnorm(n/2, -6)
        error2 <- rnorm(n/2, 6) ## shifted to create better separation :)
        index <- sample(n, n, replace = FALSE)
        error <- c(error1, error2)[index]
        y <- 5 + 2 * x1 + error
        
        dat_train <- data.frame(y = y[1:1000], x1 = x1[1:1000])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x1 = x1[(ntest + ntrain + 1):n])
        
        fit <- fit_predictionBands(dat_train$x1, dat_train$y)
        
        bands <- predict(fit, as.matrix(dat_out$x1), type = "hpd", alpha = 0.1)
        
        temp_cov <- rep(NA, length(dat_out$y))
        temp_len <- rep(NA, length(dat_out$y))
        for(j in seq(length(dat_out$y))){
            keep <- c(FALSE, unlist(bands$prediction_bands_which_belong[j]), FALSE)
            grids <- c(min(y), bands$y_grid, max(y))
            interval <- grids[keep]
            if(sum(abs(diff(keep))) > 2 & sum(abs(diff(keep))) < 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] 
            } else if(sum(abs(diff(keep))) > 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]] |
                    dat_out$y[j] > grids[index[5]] & dat_out$y[j] < grids[max(index)]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] + grids[max(index)] - grids[index[5]]
            } else{
                temp_cov[j] <- dat_out$y[j] > min(interval) & dat_out$y[j] < max(interval)
                temp_len[j] <- max(interval) - min(interval)
            }
        }
        
        cov2[i, ] <- (temp_cov) #quantile regression and upper/lower adjustments
        len2[i, ] <- (temp_len)
        
        print(i)
    }
)

mean(cov2)
sd(rowMeans(cov2)) / sqrt(N)
mean(len2)
sd(rowMeans(len2)) / sqrt(N)



## Here we look at a skewed error term, dependent on X (heteroskedastic)

cov2_uni <- matrix(NA, nrow = N, ncol = 50) 
len2_uni <- matrix(NA, nrow = N, ncol = 50) 

set.seed(6)
time_theirs_uni <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500
        x1 <- runif(n, -5, 5)
        error <- rgamma(n, shape = 1 + 2 * abs(x1), rate = 1 + 2 * abs(x1))
        y <- 5 + 2 * x1 + error
        
        dat_train <- data.frame(y = y[1:1000], x1 = x1[1:1000])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x1 = x1[(ntest + ntrain + 1):n])
        
        fit <- fit_predictionBands(dat_train$x1, dat_train$y)
        
        bands <- predict(fit, as.matrix(dat_out$x1), type = "hpd", alpha = 0.1)
        
        temp_cov <- rep(NA, length(dat_out$y))
        temp_len <- rep(NA, length(dat_out$y))
        for(j in seq(length(dat_out$y))){
            keep <- c(FALSE, unlist(bands$prediction_bands_which_belong[j]), FALSE)
            grids <- c(min(y), bands$y_grid, max(y))
            interval <- grids[keep]
            if(sum(abs(diff(keep))) > 2 & sum(abs(diff(keep))) < 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] 
            } else if(sum(abs(diff(keep))) > 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]] |
                    dat_out$y[j] > grids[index[5]] & dat_out$y[j] < grids[max(index)]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] + grids[max(index)] - grids[index[5]]
            } else{
                temp_cov[j] <- dat_out$y[j] > min(interval) & dat_out$y[j] < max(interval)
                temp_len[j] <- max(interval) - min(interval)
            }
        }
        
        cov2_uni[i, ] <- (temp_cov) #quantile regression and upper/lower adjustments
        len2_uni[i, ] <- (temp_len)
        
        print(i)
    }
)

mean(cov2_uni)
sd(rowMeans(cov2_uni)) / sqrt(N)
mean(len2_uni)
sd(rowMeans(len2_uni)) / sqrt(N)


# bowtie
cov2_bow <- matrix(NA, nrow = N, ncol = 50) 
len2_bow <- matrix(NA, nrow = N, ncol = 50) 

set.seed(6)
time_theirs_bow <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500
        x <- runif(n, -5, 5)
        error <- rnorm(n, 0, abs(x))
        y <- 5 + 2 * x + error
        
        dat_train <- data.frame(y = y[1:1000], x = x[1:1000])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x = x[(ntest + ntrain + 1):n])
        
        fit <- fit_predictionBands(dat_train$x, dat_train$y)
        
        bands <- predict(fit, as.matrix(dat_out$x), type = "hpd", alpha = 0.1)
        
        temp_cov <- rep(NA, length(dat_out$y))
        temp_len <- rep(NA, length(dat_out$y))
        for(j in seq(length(dat_out$y))){
            keep <- c(FALSE, unlist(bands$prediction_bands_which_belong[j]), FALSE)
            grids <- c(min(y), bands$y_grid, max(y))
            interval <- grids[keep]
            if(sum(abs(diff(keep))) > 2 & sum(abs(diff(keep))) < 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] 
            } else if(sum(abs(diff(keep))) > 5){
                index <- which(diff(keep) != 0) ##assuming that the first change is False to True
                temp_cov[j] <- dat_out$y[j] > grids[index[1]] & dat_out$y[j] < grids[index[2]] | 
                    dat_out$y[j] > grids[index[3]] & dat_out$y[j] < grids[index[4]] |
                    dat_out$y[j] > grids[index[5]] & dat_out$y[j] < grids[max(index)]
                temp_len[j] <- grids[index[2]] - grids[index[1]] + grids[index[4]] - grids[index[3]] + grids[max(index)] - grids[index[5]]
            } else{
                temp_cov[j] <- dat_out$y[j] > min(interval) & dat_out$y[j] < max(interval)
                temp_len[j] <- max(interval) - min(interval)
            }
        }
        
        cov2_bow[i, ] <- temp_cov #quantile regression and upper/lower adjustments
        len2_bow[i, ] <- (temp_len)
        
        print(i)
    }
)
mean(cov2_bow)
sd(rowMeans(cov2_bow)) / sqrt(N)
mean(len2_bow)
sd(rowMeans(len2_bow)) / sqrt(N)


save.image(file = "hpd-split_simulations_run.RData")

# boxplot_dat1 <- cbind(len, len2)
# colnames(boxplot_dat1) <- c("KDE-HPD", "HPD-Split")
# pdf("boxplot_comparison_bimodal.pdf", width = 8, height = 4)
# boxplot(boxplot_dat1, beside = TRUE, horizontal = TRUE, main = "Length Comparisons-Bimodal")
# dev.off()
# 
# boxplot_dat2 <- cbind(len_uni[len_uni > 1], len2_uni[len_uni > 1], len_adj)
# colnames(boxplot_dat2) <- c("KDE-HPD", "HPD-Split", "KDE-HPD-ADJ")
# pdf("boxplot_comparison_heteroskedastic.pdf", width = 8, height = 4)
# boxplot(boxplot_dat2, beside = TRUE, horizontal = TRUE, main = "Length Comparisons-Heteroskedastic")
# dev.off()
# 
# boxplot_dat3 <- cbind(len_uni_ind[len_uni_ind > 5], len2_uni_ind[len_uni_ind > 5])
# colnames(boxplot_dat3) <- c("KDE-HPD", "HPD-Split")
# pdf("boxplot_comparison_skewed.pdf", width = 8, height = 4)
# boxplot(boxplot_dat3, beside = TRUE, horizontal = TRUE, main = "Length Comparisons-Skewed")
# dev.off()
# 
# boxplot_dat4 <- cbind(len_uni_norm, len2_uni_norm)
# colnames(boxplot_dat4) <- c("KDE-HPD", "HPD-Split")
# pdf("boxplot_comparison_normal.pdf", width = 8, height = 4)
# boxplot(boxplot_dat4, beside = TRUE, horizontal = TRUE, main = "Length Comparisons-Normal")
# dev.off()
