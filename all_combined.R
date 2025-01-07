rm(list = ls())

require(hdrcde)


quant_fun <- function(dat, alpha, n_size, direction){
    o <- order(dat)
    if(all.equal(alpha, 1) == TRUE) {quantile(dat, 1)}
    else if(all.equal(alpha, 0) == TRUE) {quantile(dat, 0)}
    else if(direction != "low"){dat[o][ceiling((alpha) * (n_size + 1))]}
    else if(direction == "low"){dat[o][ceiling((alpha) * (n_size + 1)) - 1]}
}
score_fun <- function(y_hat, y_true){ 
    y_true - y_hat
}
N <- 1e3

## here we look at a standard normal error (unimodal and symmetric)

cov_uni_norm <- matrix(NA, nrow = N, ncol = 50) 
len_uni_norm <- matrix(NA, nrow = N, ncol = 50) 
set.seed(6)
time_ours_norm <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500
        x1 <- runif(n, -5, 5)
        error <- rnorm(n)
        y <- 5 + 2 * x1 + error
        
        dat_train <- data.frame(y = y[1:ntrain], x1 = x1[1:ntrain])
        dat_test <- data.frame(y = y[(ntrain + 1):(ntest + ntrain)], x1 = x1[(ntrain + 1):(ntest + ntrain)])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x1 = x1[(ntest + ntrain + 1):n])
        
        mod_train <- lm(y ~ x1, data = dat_train)
        
        pred_test <- predict(mod_train, newdata = dat_test)
        
        upper_scores <- score_fun(pred_test, dat_test$y)
        lower_scores <- -1 * upper_scores
        
        dens <- density(lower_scores, kernel = "gaussian", from = range(lower_scores)[1], to = range(lower_scores)[2], n = 500, bw = "nrd0", adjust = 500 ^(-1/3) / 500 ^(-1/5))
        hpd_result <- hdr(lower_scores, prob = 90, den = dens)
        
        cdf <- ecdf(lower_scores)
        if(length(hpd_result$hdr) > 2){
            
            lower_alpha1 <- cdf(max(hpd_result$hdr)) - cdf(hpd_result$hdr[3]) + cdf(hpd_result$hdr[2]) - 0.9
            lower_alpha2 <- cdf(hpd_result$hdr[3])
            
            upper_alpha1 <- cdf(hpd_result$hdr[2]) 
            upper_alpha2 <- cdf(max(hpd_result$hdr)) 
            
            upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ntest, "upper")
            lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ntest, "low")
            
            upper_error2 <- quant_fun(lower_scores, 1-upper_alpha2, ntest, "upper")
            lower_error2 <- quant_fun(lower_scores, 1-lower_alpha2, ntest, "low")
            
            pred_out <- predict(mod_train, newdata = dat_out)
            
            upper_out1 <- pred_out - upper_error1 
            lower_out1 <- pred_out - lower_error1 
            upper_out2 <- pred_out - upper_error2 
            lower_out2 <- pred_out - lower_error2
            cov_uni_norm[i, ] <- (dat_out$y > lower_out1 & dat_out$y < upper_out1 | dat_out$y > lower_out2 & dat_out$y < upper_out2) 
            len_uni_norm[i, ] <- (upper_out1 - lower_out1 + upper_out2 - lower_out2)
        }
        else{
            lower_alpha1 <- cdf(max(hpd_result$hdr)) - 0.9 ## so the error adds up to 1-\alpha
            
            upper_alpha1 <- cdf(max(hpd_result$hdr)) 
            
            upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ntest, "upper")
            lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ntest, "low")
            
            pred_out <- predict(mod_train, newdata = dat_out)
            
            upper_out1 <- pred_out - upper_error1 
            lower_out1 <- pred_out - lower_error1 
            cov_uni_norm[i, ] <- (dat_out$y > lower_out1 & dat_out$y < upper_out1)
            len_uni_norm[i, ] <- (upper_out1 - lower_out1)
        }
        print(i)
    }
    
)

# plot(y ~ x1, ylab = "Y", xlab = "X", main = "", col = "gray")
# fill_points <- seq(from = -5, to = 5, length.out = 1000)
# graph_preds <- predict(mod_train, newdata = data.frame(x1 = fill_points))
# y_plot_low1 <- graph_preds - lower_error1
# y_plot_high1 <- graph_preds - upper_error1
# polygon(x = c(-5, -5, 5, 5), y = c(y_plot_low1[1], y_plot_high1[1], y_plot_high1[1000], y_plot_low1[1000]), col = rgb(0, 0, 0, 0.5))
mean(cov_uni_norm)
sd(rowMeans(cov_uni_norm)) / sqrt(N)
mean(len_uni_norm)
sd(rowMeans(len_uni_norm)) / sqrt(N)

## here we look at a skewed error independent of X (unimodal and skewed)
cov_uni_ind <- matrix(NA, nrow = N, ncol = 50) 
len_uni_ind <- matrix(NA, nrow = N, ncol = 50) 
set.seed(6)
time_ours_ind <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500

        x1 <- runif(n, -5, 5)
        error <- rgamma(n, shape = 7.5, rate = 1 )
        y <- 5 + 2 * x1 + error
        
        dat_train <- data.frame(y = y[1:ntrain], x1 = x1[1:ntrain])
        dat_test <- data.frame(y = y[(ntrain + 1):(ntest + ntrain)], x1 = x1[(ntrain + 1):(ntest + ntrain)])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x1 = x1[(ntest + ntrain + 1):n])
        
        mod_train <- lm(y ~ x1, data = dat_train)
        
        pred_test <- predict(mod_train, newdata = dat_test)
        
        upper_scores <- score_fun(pred_test, dat_test$y)
        lower_scores <- -1 * upper_scores


        dens <- density(lower_scores, kernel = "gaussian", from = range(lower_scores)[1], to = range(lower_scores)[2], n = 500, bw = "nrd0", adjust = 500 ^(-1/3) / 500 ^(-1/5))
        hpd_result <- hdr(lower_scores, prob = 90, den = dens)
        cdf <- ecdf(lower_scores)
        if(length(hpd_result$hdr) > 2){
            
            lower_alpha1 <- cdf(max(hpd_result$hdr)) - cdf(hpd_result$hdr[3]) + cdf(hpd_result$hdr[2]) - 0.9
            lower_alpha2 <- cdf(hpd_result$hdr[3])
            
            upper_alpha1 <- cdf(hpd_result$hdr[2]) 
            upper_alpha2 <- cdf(max(hpd_result$hdr)) 
            
            upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ntest, "upper")
            lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ntest, "low")

            upper_error2 <- quant_fun(lower_scores, 1-upper_alpha2, ntest, "upper")
            lower_error2 <- quant_fun(lower_scores, 1-lower_alpha2, ntest, "low")

            pred_out <- predict(mod_train, newdata = dat_out)
            
            upper_out1 <- pred_out - upper_error1 
            lower_out1 <- pred_out - lower_error1 
            upper_out2 <- pred_out - upper_error2 
            lower_out2 <- pred_out - lower_error2
            cov_uni_ind[i, ] <- (dat_out$y > lower_out1 & dat_out$y < upper_out1 | dat_out$y > lower_out2 & dat_out$y < upper_out2) 
            len_uni_ind[i, ] <- (upper_out1 - lower_out1 + upper_out2 - lower_out2)
        }
        else{
            lower_alpha1 <- cdf(max(hpd_result$hdr)) - 0.9
            
            upper_alpha1 <- cdf(max(hpd_result$hdr)) 
            
            upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ntest, "upper")
            lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ntest, "low")

            pred_out <- predict(mod_train, newdata = dat_out)
            
            upper_out1 <- pred_out - upper_error1 
            lower_out1 <- pred_out - lower_error1 
            cov_uni_ind[i, ] <- (dat_out$y > lower_out1 & dat_out$y < upper_out1) 
            len_uni_ind[i, ] <- (upper_out1 - lower_out1)
        }
        print(i)
    }
)

# plot(y ~ x1, ylab = "Y", xlab = "X", main = "", col = "gray")
# fill_points <- seq(from = -5, to = 5, length.out = 1000)
# graph_preds <- predict(mod_train, newdata = data.frame(x1 = fill_points))
# y_plot_low1 <- graph_preds - lower_error1
# y_plot_high1 <- graph_preds - upper_error1
# polygon(x = c(-5, -5, 5, 5), y = c(y_plot_low1[1], y_plot_high1[1], y_plot_high1[1000], y_plot_low1[1000]), col = rgb(0, 0, 0, 0.5))


mean(cov_uni_ind)
sd(rowMeans(cov_uni_ind)) / sqrt(N)
mean(len_uni_ind)
sd(rowMeans(len_uni_ind)) / sqrt(N)

## here we look at a bimodal error term, independent of X (bimodal)

cov <- matrix(NA, nrow = N, ncol = 50) 
len <- matrix(NA, nrow = N, ncol = 50) 
set.seed(6)
time_ours <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500
        
        x1 <- runif(n, -5, 5)
        index <- sample(n, n, replace = FALSE)

        error1 <- rnorm(n/2, -6)
        error2 <- rnorm(n/2, 6)
        error <- c(error1, error2)[index]
        y <- 5 + 2 * x1 + error
        
        dat_train <- data.frame(y = y[1:ntrain], x1 = x1[1:ntrain])
        dat_test <- data.frame(y = y[(ntrain + 1):(ntest + ntrain)], x1 = x1[(ntrain + 1):(ntest + ntrain)])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x1 = x1[(ntest + ntrain + 1):n])
        
        mod_train <- lm(y ~ x1, data = dat_train)
        
        pred_test <- predict(mod_train, newdata = dat_test)
        
        upper_scores <- score_fun(pred_test, dat_test$y)
        lower_scores <- -1 * upper_scores
        
        
        
        dens <- density(lower_scores, kernel = "gaussian", from = range(lower_scores)[1], to = range(lower_scores)[2], n = 500, bw = "nrd0", adjust = 500 ^(-1/3) / 500 ^(-1/5))
        hpd_result <- hdr(lower_scores, prob = 90, den = dens)
        
        
        cdf <- ecdf(lower_scores)

        lower_alpha1 <- cdf(max(hpd_result$hdr)) - cdf(hpd_result$hdr[3]) + cdf(hpd_result$hdr[2]) - 0.9
        lower_alpha2 <- cdf(hpd_result$hdr[3])
        
        upper_alpha1 <- cdf(hpd_result$hdr[2]) 
        upper_alpha2 <- cdf(max(hpd_result$hdr)) 
        
        upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ntest, "upper")
        lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ntest, "low")
        
        upper_error2 <- quant_fun(lower_scores, 1-upper_alpha2, ntest, "upper")
        lower_error2 <- quant_fun(lower_scores, 1-lower_alpha2, ntest, "low")
        
        
        pred_out <- predict(mod_train, newdata = dat_out)
        
        upper_out1 <- pred_out - upper_error1
        lower_out1 <- pred_out - lower_error1
        upper_out2 <- pred_out - upper_error2
        lower_out2 <- pred_out - lower_error2
        cov[i, ] <- (dat_out$y > lower_out1 & dat_out$y < upper_out1 | dat_out$y > lower_out2 & dat_out$y < upper_out2) 
        if(upper_out1[1] > lower_out2[1]){len[i, ] <- (upper_out2 - lower_out1)}
        else {len[i, ] <- (upper_out1 - lower_out1 + upper_out2 - lower_out2)}
        print(i)
    }
)
# plot(y[1:5000] ~ x1[1:5000], ylab = "Y", xlab = "X", main = "", col = "gray")
# fill_points <- seq(from = -5, to = 5, length.out = 1000)
# graph_preds <- predict(mod_train, newdata = data.frame(x1 = fill_points))
# y_plot_low1 <- graph_preds - lower_error1
# y_plot_high1 <- graph_preds - upper_error1
# y_plot_low2 <- graph_preds - lower_error2
# y_plot_high2 <- graph_preds - upper_error2
# polygon(x = c(-5, -5, 5, 5), y = c(y_plot_low1[1], y_plot_high1[1], y_plot_high1[1000], y_plot_low1[1000]), col = rgb(0, 0, 0, 0.5))
# polygon(x = c(-5,-5 , 5, 5), y = c(y_plot_low2[1], y_plot_high2[1], y_plot_high2[1000], y_plot_low2[1000]), col = rgb(0, 0, 0, 0.5))


mean(cov)
sd(rowMeans(cov)) / sqrt(N)
mean(len)
sd(rowMeans(len)) / sqrt(N)


## (heteroskedastic, no adj)

cov_adj <- matrix(NA, nrow = N, ncol = 50) 
len_adj <- matrix(NA, nrow = N, ncol = 50) 
set.seed(6)
time_ours_adj <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntest <- 500
        ntrain <- 500
        x1 <- runif(n, -5, 5)
        error <- rgamma(n, shape = 1 + 2 * abs(x1), rate = 1 + 2 * abs(x1))
        y <- 5 + 2 * x1 + error
        
        dat_train <- data.frame(y = y[1:ntrain], x1 = x1[1:ntrain])
        dat_test <- data.frame(y = y[(ntrain + 1):(ntest + ntrain)], x1 = x1[(ntrain + 1):(ntest + ntrain)])
        dat_out <- data.frame(y = y[(ntest + ntrain + 1):n], x1 = x1[(ntest + ntrain + 1):n])
        
        mod_train <- lm(y ~ x1, data = dat_train)
        
        pred_test <- predict(mod_train, newdata = dat_test)
        
        upper_scores <- score_fun(pred_test, dat_test$y)
        lower_scores <- -1 * upper_scores
        
        dens <- density(lower_scores, kernel = "gaussian", from = range(lower_scores)[1], to = range(lower_scores)[2], n = 500, bw = "nrd0", adjust = 500 ^(-1/3) / 500 ^(-1/5))
        hpd_result <- hdr(lower_scores, prob = 90, den = dens)
        
        cdf <- ecdf(lower_scores)
        
        if(length(hpd_result$hdr) > 3){
            
            lower_alpha1 <- cdf(max(hpd_result$hdr)) - cdf(hpd_result$hdr[3]) + cdf(hpd_result$hdr[2]) - 0.9
            lower_alpha2 <- cdf(hpd_result$hdr[3])
            
            upper_alpha1 <- cdf(hpd_result$hdr[2]) 
            upper_alpha2 <- cdf(max(hpd_result$hdr)) 
            
            upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ntest, "upper")
            lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ntest, "low")
            
            upper_error2 <- quant_fun(lower_scores, 1-upper_alpha2, ntest, "upper")
            lower_error2 <- quant_fun(lower_scores, 1-lower_alpha2, ntest, "low")
            
            pred_out <- predict(mod_train, newdata = dat_out)

            upper_out1 <- pred_out - upper_error1 
            lower_out1 <- pred_out - lower_error1 
            upper_out2 <- pred_out - upper_error2 
            lower_out2 <- pred_out - lower_error2
            cov_adj[i, ] <- (dat_out$y > lower_out1 & dat_out$y < upper_out1 | dat_out$y > lower_out2 & dat_out$y < upper_out2) 
            len_adj[i, ] <- (upper_out1 - lower_out1 + upper_out2 - lower_out2)
        }
        else{
            lower_alpha1 <- cdf(max(hpd_result$hdr)) - 0.9
            
            upper_alpha1 <- cdf(max(hpd_result$hdr)) 
            
            upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ntest, "upper")
            lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ntest, "low")
            
            pred_out <- predict(mod_train, newdata = dat_out)

            upper_out1 <- pred_out - upper_error1 
            lower_out1 <- pred_out - lower_error1 
            cov_adj[i, ] <- (dat_out$y > lower_out1 & dat_out$y < upper_out1) 
            len_adj[i, ] <- (upper_out1 - lower_out1)
        }
        print(i)
    }
)
# plot(y ~ x1, ylab = "Y", xlab = "X", main = "", col = "gray")
# fill_points <- seq(from = -5, to = 5, length.out = 1000)
# graph_preds <- predict(mod_train, newdata = data.frame(x1 = fill_points))
# y_plot_low1 <- graph_preds - lower_error1
# y_plot_high1 <- graph_preds - upper_error1
# polygon(x = c(-5, -5, 5, 5), y = c(y_plot_low1[1], y_plot_high1[1], y_plot_high1[1000], y_plot_low1[1000]), col = rgb(0, 0, 0, 0.5))

mean(cov_adj)
sd(rowMeans(cov_adj)) / sqrt(N)
mean(len_adj)
sd(rowMeans(len_adj)) / sqrt(N)


# bowtie
cov_bow <- matrix(NA, nrow = N, ncol = 50) 
len_bow <- matrix(NA, nrow = N, ncol = 50) 
set.seed(6)
time_ours_bow <- system.time(
    for(i in seq(N)){
        n <- 1050
        ntrain <- 250
        nhetero <- 250
        ncal <- 500
        x <- runif(n, -5, 5)
        error <- rnorm(n, 0, abs(x))
        y <- 5 + 2 * x + error
        
        dat_train <- data.frame(y = y[1:ntrain], x = x[1:ntrain])
        dat_hetero <- data.frame(y = y[251:500], x = x[251:500])
        dat_cal <- data.frame(y = y[(500 + 1):1000], x = x[501:1000])
        dat_out <- data.frame(y = y[1001:n], x = x[1001:n])
        
        mod_train <- lm(y ~ x, data = dat_train)
        
        hetero_test <- predict(mod_train, newdata = dat_hetero)
        sigma_test <- abs(hetero_test - dat_hetero$y)
        dat_hetero <- cbind(dat_hetero, sigma_test)
        
        mod_sig <- quantregRanger::quantregRanger(sigma_test ~ x, data = dat_hetero)
        
        sig_cal <- predict(mod_sig, data = dat_cal, quantiles = 0.9)
        sig_out <- predict(mod_sig, data = dat_out, quantiles = 0.9)
        
        pred_cal <- predict(mod_train, newdata = dat_cal)
        
        upper_scores <- score_fun(pred_cal, dat_cal$y) / abs(sig_cal)
        lower_scores <- -1 * upper_scores
        
        
        dens <- density(lower_scores, kernel = "gaussian", from = range(lower_scores)[1], to = range(lower_scores)[2], n = 500, adjust = 500 ^(-1/3) / 500 ^(-1/5))
        hpd_result <- hdr(lower_scores, prob = 90, den = dens)
        cdf <- ecdf(lower_scores)
        if(length(hpd_result$hdr) > 2){
            lower_alpha1 <- cdf(max(hpd_result$hdr)) - cdf(hpd_result$hdr[3]) + cdf(hpd_result$hdr[2]) - 0.9
            lower_alpha2 <- cdf(hpd_result$hdr[3])
            
            upper_alpha1 <- cdf(hpd_result$hdr[2]) 
            upper_alpha2 <- cdf(max(hpd_result$hdr)) 
            
            upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ncal, "upper")
            lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ncal, "low")
            
            upper_error2 <- quant_fun(lower_scores, 1-upper_alpha2, ncal, "upper")
            lower_error2 <- quant_fun(lower_scores, 1-lower_alpha2, ncal, "low")
            
            pred_out <- predict(mod_train, newdata = dat_out)
            
            upper_out1 <- pred_out - upper_error1 * abs(sig_out)
            lower_out1 <- pred_out - lower_error1 * abs(sig_out)
            upper_out2 <- pred_out - upper_error2 * abs(sig_out)
            lower_out2 <- pred_out - lower_error2 * abs(sig_out)
            cov_bow[i, ] <- dat_out$y > lower_out1 & dat_out$y < upper_out1 | dat_out$y > lower_out2 & dat_out$y < upper_out2
            len_bow[i, ] <- (upper_out1 - lower_out1 + upper_out2 - lower_out2)
            print("hi")
        }
        else{
            lower_alpha1 <- cdf(max(hpd_result$hdr)) - 0.9
            
            upper_alpha1 <- cdf(max(hpd_result$hdr)) 
            
            upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, ncal, "upper")
            lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, ncal, "low")
            
            pred_out <- predict(mod_train, newdata = dat_out)
            
            upper_out1 <- (pred_out - upper_error1 * abs(sig_out))
            lower_out1 <- (pred_out - lower_error1 * abs(sig_out))
            cov_bow[i, ] <- dat_out$y > lower_out1 & dat_out$y < upper_out1
            len_bow[i, ] <- upper_out1 - lower_out1
        }
        print(i)
        
    }
)

mean(cov_bow)
sd(rowMeans(cov_bow)) / sqrt(N)
mean(len_bow)
sd(rowMeans(len_bow)) / sqrt(N)


# plot(y ~ x, ylab = "Y", xlab = "X", main = "", col = "gray")
# fill_points <- seq(from = -5, to = 5, length.out = 1000)
# graph_preds <- predict(mod_train, newdata = data.frame(x = fill_points))
# sig_graph <- predict(mod_sig, data = data.frame(x = fill_points), quantiles = 0.9)
# 
# y_plot_low1 <- graph_preds - lower_error1 * abs(sig_graph)
# y_plot_high1 <- graph_preds - upper_error1 * abs(sig_graph)
# y_plot_combined <- c(y_plot_low1, y_plot_high1)[order(c(seq_along(y_plot_low1), seq_along(y_plot_high1)))]
# #polygon(x = rep(fill_points, 2), y = c(y_plot_low1, y_plot_high1), col = rgb(0, 0, 0, 0.5))
# polygon(x = c(rbind(fill_points, fill_points)), y = y_plot_combined, col = rgb(0, 0, 0, 0.5))

# x_out <- seq(-5.07, 5.07, length.out = 50)
# graph_preds <- predict(mod_train, newdata = data.frame(x = x_out))
# sig_graph <- predict(mod_sig, data = data.frame(x = x_out), quantiles = 0.9)
# y_plot_low1 <- as.vector(graph_preds - lower_error1 * sig_graph)
# y_plot_high1 <- as.vector(graph_preds - upper_error1 * sig_graph)
# plot_dat <- data.frame(x = x, y = y)
# plot_dat2 <- data.frame(ylow1 = y_plot_low1, yhigh1 = y_plot_high1, x = x_out)
# ggplot(data = plot_dat2, aes(x = x)) + 
#     geom_ribbon(data = plot_dat2, aes(ymin = ylow1, ymax = yhigh1), alpha = 1, fill = "darkgrey") + 
#     geom_point(data = plot_dat, aes(x = x, y = y), shape = 21, alpha = 0.5, fill = NA) +
#     theme_classic() 


#save.image(file = "kde-hpd_simulations_run.RData")

